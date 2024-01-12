! =====================================================================
	subroutine QQ2xx_SetCounterterms()
	implicit none

#include "QQ2xx_Model.h"
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer i,j,k,l,m,type,gen
        integer typex,typey, sfmeigenin, sfmeigenout, sfgaugeeigenin, sfgaugeeigenout
	integer ttf, genf,ttsf,ttsfa,ttsfb,gensf,gensfa,gensfb, eigenchi, eigensf
	double precision CL1, CL2, Krondelta
        double precision CLR11,CLR22
        double complex dCLnew

	external QQ2xx_CalcPropSquark, QQ2xx_CalcPropQuark

! ===================================
! ===== Coupling counterterms ======= 
! ===================================

! =====================================================================================================		
! ================================== Vertex Squark - Squark - Vector ===================================
! =====================================================================================================	
!see VSfSf

        do sfmeigenout = 1, 2
          do eigensf = 1, 2
            do ttsf = 3, 4
              do gensf = 3, 3

              CL1 = - g2h1*I3f_(ttsf)/CW
              CL2 = g2h1*SW**2*ef_(ttsf)/CW

        
              dVSfSf(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf,gensf)) +
     &        (CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf,gensf) + CL2*Krondelta(eigensf,1))*dconjg(dZSf(1,sfmeigenout,ttsf,gensf))/2. +
     &        (CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf,gensf) + CL2*Krondelta(eigensf,2))*dconjg(dZSf(2,sfmeigenout,ttsf,gensf))/2. +
     &        (CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + CL2*Krondelta(1,sfmeigenout))*dZSf(1,eigensf,ttsf,gensf)/2. +
     &        (CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + CL2*Krondelta(2,sfmeigenout))*dZSf(2,eigensf,ttsf,gensf)/2.

              dVSfSft(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = dVSfSf(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) 
              dVSfSfu(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = dVSfSf(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) 

     		   CLR11= -g2h1/CW*(I3f_(ttsf)-SW**2*ef_(ttsf))
     		   CLR22= -g2h1/CW*(-SW**2*ef_(ttsf))

     		  dCLnew = 0.5d0*VSfSf(3,1,eigensf,ttsf,ttsf,gensf,gensf)*dconjg(dZSf(1,sfmeigenout,ttsf,gensf))
     &		          +0.5d0*VSfSf(3,2,eigensf,ttsf,ttsf,gensf,gensf)*dconjg(dZSf(2,sfmeigenout,ttsf,gensf))
     &   		      +0.5d0*VSfSf(3,sfmeigenout,1,ttsf,ttsf,gensf,gensf)*dZSf(1,eigensf,ttsf,gensf)
     &   		      +0.5d0*VSfSf(3,sfmeigenout,2,ttsf,ttsf,gensf,gensf)*dZSf(2,eigensf,ttsf,gensf)
     &   		      +dUSf(sfmeigenout,1,ttsf,gensf)*CLR11*USf(eigensf,1,ttsf,gensf)
     &   		      +dUSf(sfmeigenout,2,ttsf,gensf)*CLR22*USf(eigensf,2,ttsf,gensf)
     &   		      +USf(sfmeigenout,1,ttsf,gensf)*CLR11*dUSf(eigensf,1,ttsf,gensf)
     &   		      +USf(sfmeigenout,2,ttsf,gensf)*CLR22*dUSf(eigensf,2,ttsf,gensf)

c     		   print*,'dZ11=',dZSf(1,1,ttsf,gensf)
c     		   print*,'dZ22=',dZSf(2,2,ttsf,gensf)
c     		   print*,'dZoffsym=',dZSf(1,2,ttsf,gensf)+dZSf(2,1,ttsf,gensf)
c     		   print*,'dZoffas=', dZSf(1,2,ttsf,gensf)-dZSf(2,1,ttsf,gensf)
c     		   print*,'dUSf=',dUSf(eigensf,1,ttsf,gensf)
c     		   print*,'dCLQ=',dVSfSf(3,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf)
c     		   print*,'dCLK=',dCLnew
c     
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
              CL2 = -ef_(ttsf)*EL
        
              dVSfSf(2,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = CL2*(dconjg(dZSf(eigensf,sfmeigenout,ttsf,gensf)) + dZSf(sfmeigenout,eigensf,ttsf,gensf))/2.

              dVSfSft(2,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = dVSfSf(2,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) 
              dVSfSfu(2,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) = dVSfSf(2,sfmeigenout,eigensf,ttsf,ttsf,gensf,gensf) 
             
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

              dVSfSf(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf+1,gensf)) +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf+1,gensf)*dconjg(dZSf(1,sfmeigenout,ttsf+1,gensf))/2. +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf+1,gensf)*dconjg(dZSf(2,sfmeigenout,ttsf+1,gensf))/2. +
     &        CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf)*dZSf(1,eigensf,ttsf,gensf)/2. +
     &        CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf)*dZSf(2,eigensf,ttsf,gensf)/2.

               if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then  !u-channel only needs to be switched of for +- final-states (e.g. W+H- etc.)
                dVSfSft(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf)
                dVSfSfu(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf) = (0d0,0d0)    
               else
                dVSfSft(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf)  
                dVSfSfu(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf+1,ttsf,gensf,gensf)
               endif
        
! --- incoming down-type sfermion ---

              else

              dVSfSf(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf-1,gensf)) +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf-1,gensf)*dconjg(dZSf(1,sfmeigenout,ttsf-1,gensf))/2. +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf-1,gensf)*dconjg(dZSf(2,sfmeigenout,ttsf-1,gensf))/2. +
     &        CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf)*dconjg(dZSf(1,eigensf,ttsf,gensf))/2. +
     &        CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf)*dconjg(dZSf(2,eigensf,ttsf,gensf))/2.

               if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then  !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
                dVSfSft(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf) = (0d0,0d0)    
                dVSfSfu(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf)
               else
                dVSfSft(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf)    
                dVSfSfu(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf) = dVSfSf(4,sfmeigenout,eigensf,ttsf-1,ttsf,gensf,gensf)
               endif
        
              endif
        
              enddo
            enddo
          enddo
        enddo
        
! =====================================================================================================		
! ================================== Vertex Squark - Squark - Higgs ===================================
! =====================================================================================================	
!see HSfSf
        
         do ttsfa = 1, 4
           do ttsfb = 1, 4
             do gensfa = 1, 3
               do gensfb = 1, 3
                 do sfmeigenin = 1, 2
                   do sfmeigenout = 1,2
                     do typex = 1, 6

           dHSfSf(typex,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensfa,gensfb) = (0d0,0d0)

                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         enddo

     	do ttsf = 1, 4
           do gensf = 1, 3
           
! --- incoming up-type sfermion --- 

	if ((ttsf.eq.1).or.(ttsf.eq.3)) then 
! h
	dGLR(1,1,1,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,2,2,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,2,1,ttsf,gensf) = (-hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) - (hf_(ttsf,gensf)/sqrt2*MUE*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,1,2,ttsf,gensf) = dGLR(1,2,1,ttsf,gensf)

! H
	dGLR(2,1,1,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,2,2,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,2,1,ttsf,gensf) = (-hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)/sqrt2*MUE*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,1,2,ttsf,gensf) =  dGLR(2,2,1,ttsf,gensf)
	
! H3
	dGLR(3,1,1,ttsf,gensf) = 0d0
	dGLR(3,2,2,ttsf,gensf) = 0d0
	dGLR(3,2,1,ttsf,gensf) = (cI*hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (cI*hf_(ttsf,gensf)/sqrt2*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(3,1,2,ttsf,gensf) = - dGLR(3,2,1,ttsf,gensf)
        
! G0
       	dGLR(4,1,1,ttsf,gensf) = 0d0
	dGLR(4,2,2,ttsf,gensf) = 0d0
	dGLR(4,2,1,ttsf,gensf) = (cI*hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) - (cI*hf_(ttsf,gensf)/sqrt2*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(4,1,2,ttsf,gensf) = - dGLR(4,2,1,ttsf,gensf)

! H+
	dGLR(5,1,1,ttsf,gensf) = (hf_(ttsf+1,gensf)*Mf(ttsf+1,gensf)*SB)*(2d0*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf)) + (hf_(ttsf,gensf)*Mf(ttsf,gensf)*CB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,2,2,ttsf,gensf) = ((hf_(ttsf,gensf)*Mf(ttsf+1,gensf)*CB) + hf_(ttsf+1,gensf)*Mf(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))
	dGLR(5,2,1,ttsf,gensf) = (hf_(ttsf,gensf)*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,1,2,ttsf,gensf) = (hf_(ttsf+1,gensf)*Af(ttsf+1,gensf)*SB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf) + dAf(ttsf+1,gensf)/Af(ttsf+1,gensf)) + (hf_(ttsf+1,gensf)*MUE*CB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))

! G+
	dGLR(6,1,1,ttsf,gensf) = -(hf_(ttsf+1,gensf)*Mf(ttsf+1,gensf)*CB)*(2d0*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf)) + (hf_(ttsf,gensf)*Mf(ttsf,gensf)*SB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,2,ttsf,gensf) = ((hf_(ttsf,gensf)*Mf(ttsf+1,gensf)*SB) - hf_(ttsf+1,gensf)*Mf(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))
	dGLR(6,2,1,ttsf,gensf) = (hf_(ttsf,gensf)*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) - (hf_(ttsf,gensf)*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,1,2,ttsf,gensf) = -(hf_(ttsf+1,gensf)*Af(ttsf+1,gensf)*CB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf) + dAf(ttsf+1,gensf)/Af(ttsf+1,gensf)) + (hf_(ttsf+1,gensf)*MUE*SB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))

! --- incoming down-type sfermion --- 

	else  
! h
	dGLR(1,1,1,ttsf,gensf) = (sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,2,2,ttsf,gensf) = (sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,2,1,ttsf,gensf) = (hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)/sqrt2*MUE*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(1,1,2,ttsf,gensf) = dGLR(1,2,1,ttsf,gensf)
! H
	dGLR(2,1,1,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,2,2,ttsf,gensf) = (-sqrt2*hf_(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,2,1,ttsf,gensf) = (-hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)/sqrt2*MUE*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(2,1,2,ttsf,gensf) = dGLR(2,2,1,ttsf,gensf)
! H3
	dGLR(3,1,1,ttsf,gensf) = 0d0
	dGLR(3,2,2,ttsf,gensf) = 0d0
	dGLR(3,2,1,ttsf,gensf) = (cI*hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (cI*hf_(ttsf,gensf)/sqrt2*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(3,1,2,ttsf,gensf) = -dGLR(3,2,1,ttsf,gensf)

! G0
        dGLR(4,1,1,ttsf,gensf) = 0d0
	dGLR(4,2,2,ttsf,gensf) = 0d0
	dGLR(4,2,1,ttsf,gensf) = -(cI*hf_(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (cI*hf_(ttsf,gensf)/sqrt2*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(4,1,2,ttsf,gensf) = -dGLR(4,2,1,ttsf,gensf)

! H-
	dGLR(5,1,1,ttsf,gensf) = (hf_(ttsf,gensf)*Mf(ttsf,gensf)*SB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf)) + (hf_(ttsf-1,gensf)*Mf(ttsf-1,gensf)*CB)*(2d0*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))
	dGLR(5,2,2,ttsf,gensf) = (hf_(ttsf-1,gensf)*Mf(ttsf,gensf)*CB + hf_(ttsf,gensf)*Mf(ttsf-1,gensf)*SB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dMf(ttsf,gensf)/Mf(ttsf,gensf))	
 	dGLR(5,2,1,ttsf,gensf) = (hf_(ttsf,gensf)*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,1,2,ttsf,gensf) = (hf_(ttsf-1,gensf)*Af(ttsf-1,gensf)*CB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dAf(ttsf-1,gensf)/Af(ttsf-1,gensf)) + (hf_(ttsf-1,gensf)*MUE*SB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))

! G-
	dGLR(6,1,1,ttsf,gensf) = -(hf_(ttsf,gensf)*Mf(ttsf,gensf)*CB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf)) + (hf_(ttsf-1,gensf)*Mf(ttsf-1,gensf)*SB)*(2d0*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))
	dGLR(6,2,2,ttsf,gensf) = (hf_(ttsf-1,gensf)*Mf(ttsf,gensf)*SB - hf_(ttsf,gensf)*Mf(ttsf-1,gensf)*CB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,1,ttsf,gensf) = -(hf_(ttsf,gensf)*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf_(ttsf,gensf)*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,1,2,ttsf,gensf) = (hf_(ttsf-1,gensf)*Af(ttsf-1,gensf)*SB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dAf(ttsf-1,gensf)/Af(ttsf-1,gensf)) - (hf_(ttsf-1,gensf)*MUE*CB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf)) 

	endif
	
		

! ------------------------------ Vertex Squark - Squark - neutral Higgs ------------------------------

	do typex=1,4
	   do sfmeigenin=1,2
	      do sfmeigenout=1,2
	         do sfgaugeeigenin=1,2
	            do sfgaugeeigenout=1,2

    	 dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) +
     &   dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*dGLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf,gensf)	               
   	            enddo

		dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) + 
     &          1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf,gensf))*HSfSf(typex,sfmeigenin,sfgaugeeigenin,ttsf,ttsf,gensf,gensf) +               
     &          dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HSfSf(typex,sfgaugeeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf))
	         enddo

        dHSfSft(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf)      
        dHSfSfu(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf,gensf) 

 	      enddo
	   enddo
	enddo



! ------------------------------ Vertex Squark - Squark - charged Higgs ------------------------------
	do typex=5,6
          do sfmeigenin=1,2
	    do sfmeigenout=1,2
! --- incoming up-type sfermion ---	
	      if ((ttsf.eq.1).or.(ttsf.eq.3)) then
	
	      do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
	            dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) +
     &              dGLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                 enddo
		    dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf+1,gensf))*HSfSf(typex,sfmeigenin,sfgaugeeigenin,ttsf,ttsf+1,gensf,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HSfSf(typex,sfgaugeeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf))
	      enddo

                if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                          !u-channel only needs to be switched of for +- final-states (e.g. W+H- etc.)
             
                  dHSfSft(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf)      !otherwise use usual couplings and switch off u-channel in amplitude wire if-statements
                  dHSfSfu(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = (0d0,0d0)
              
                else
    
                  dHSfSft(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf)
                  dHSfSfu(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf,gensf)

                endif



! --- incoming down-type sfermion ---
                    
              else

              do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
                    dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) +
     &              dGLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              GLR(typex,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                 enddo
		    dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf-1,gensf))*HSfSf(typex,sfmeigenin,sfgaugeeigenin,ttsf,ttsf-1,gensf,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HSfSf(typex,sfgaugeeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf))
 		
	      enddo

                  if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                        !u-channel only needs to be switched of for +- final-states (e.g. W+H- etc.)
                  
                    dHSfSft(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = (0d0,0d0)
                    dHSfSfu(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf)

                  else

                    dHSfSft(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf)
                    dHSfSfu(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf) = dHSfSf(typex,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf,gensf)
        
                  endif


             endif

            enddo
	  enddo
	enddo
	
	enddo
	enddo

!======================================================================================================================
!======================================  Vertex sFermion-sFermion-Higgs-Higgs =========================================
!======================================================================================================================
!see HHSfSf


        do typex=1,6                 !no flavor-mixing, therefore sfgen2 not needed
          do typey=1,6
            do sfmeigenin=1,2
              do sfmeigenout= 1,2
                do ttsfa=1,4
                  do gensf=1,3

                   dKLR(typex,typey,sfmeigenin,sfmeigenout,ttsfa,gensf) = (0d0,0d0)

                   do ttsfb=1,4

                   dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = (0d0,0d0)

                   enddo

                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do ttsf= 1, 4                    
           do gensf=1, 3

	if ((ttsf.eq.1).or.(ttsf.eq.3)) then
       
!h,h

	dKLR(1,1,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CA2                     
	dKLR(1,1,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CA2               

!H,H
 
        dKLR(2,2,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SA2   
	dKLR(2,2,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SA2   

!P,P

        dKLR(3,3,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CB2 
	dKLR(3,3,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CB2   

!G0,G0

        dKLR(4,4,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SB2 
	dKLR(4,4,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SB2  
             
!h,H

	dKLR(1,2,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2A/2d0                  !LL-Kopplung
	dKLR(1,2,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2A/2d0                  !RR-Kopplung
	dKLR(2,1,1,1,ttsf,gensf) = dKLR(1,2,1,1,ttsf,gensf)                                                              
	dKLR(2,1,2,2,ttsf,gensf) = dKLR(1,2,2,2,ttsf,gensf)                                                                          

!G0,P
      
	dKLR(3,4,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2B/2d0
	dKLR(3,4,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2B/2d0 
	dKLR(4,3,1,1,ttsf,gensf) = dKLR(3,4,1,1,ttsf,gensf)
	dKLR(4,3,2,2,ttsf,gensf) = dKLR(3,4,2,2,ttsf,gensf)

!Hp-Hm

	dKLR(5,5,1,1,ttsf,gensf) = - 2*EL**2/(MZ2*S2W**2)*Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))*TB2
	dKLR(5,5,2,2,ttsf,gensf) = - 2*EL**2*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MZ2*S2W**2*TB2)

!Gp-Gm

	dKLR(6,6,1,1,ttsf,gensf) = - 2*EL**2/(MZ2*S2W**2)*Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))
	dKLR(6,6,2,2,ttsf,gensf) = - 2*EL**2*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MZ2*S2W**2)

!Gp,Hm (Gm,Hp)

	dKLR(6,5,1,1,ttsf,gensf) = 2*EL**2/(MZ2*S2W**2)*Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))*TB
	dKLR(6,5,2,2,ttsf,gensf) =- EL**2/(2d0*SW2*MW2*TB)*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dKLR(5,6,1,1,ttsf,gensf) = dKLR(6,5,1,1,ttsf,gensf)
	dKLR(5,6,2,2,ttsf,gensf) = dKLR(6,5,2,2,ttsf,gensf)

!neutral Higgs,Hp (Hm)
!error in the HiggsHunter p.408, necessary replacements->see Couplings.F

	dKLR(1,5,1,1,ttsf,gensf) =-EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*CA*CB/SB2-Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*(-SA)*SB/CB2) 
	dKLR(1,5,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2*SBA
	dKLR(5,1,1,1,ttsf,gensf) = dKLR(1,5,1,1,ttsf,gensf)
	dKLR(5,1,2,2,ttsf,gensf) = dKLR(1,5,2,2,ttsf,gensf)

	dKLR(2,5,1,1,ttsf,gensf) =-EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*SA*CB/SB2+Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*(-CA)*SB/CB2)
	dKLR(2,5,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2*CBA
	dKLR(5,2,1,1,ttsf,gensf) = dKLR(2,5,1,1,ttsf,gensf)
	dKLR(5,2,2,2,ttsf,gensf) = dKLR(2,5,2,2,ttsf,gensf)

	dKLR(3,5,1,1,ttsf,gensf) =-cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MW2*TB2)+Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*TB2) 
	dKLR(3,5,2,2,ttsf,gensf) = (0d0,0d0)
	dKLR(5,3,1,1,ttsf,gensf) = dKLR(3,5,1,1,ttsf,gensf)
	dKLR(5,3,2,2,ttsf,gensf) = dKLR(3,5,2,2,ttsf,gensf)

	dKLR(4,5,1,1,ttsf,gensf) =-cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MW2*TB)-Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*TB) 
	dKLR(4,5,2,2,ttsf,gensf) =-cI*hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2
	dKLR(5,4,1,1,ttsf,gensf) = dKLR(4,5,1,1,ttsf,gensf)
	dKLR(5,4,2,2,ttsf,gensf) = dKLR(4,5,2,2,ttsf,gensf)

!neutral Higgs,Gp (Gm)
!error in the HiggsHunter p.408, necessary replacements->see Couplings.F

	dKLR(1,6,1,1,ttsf,gensf) = EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*(-CA)/SB+Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*SA/CB) 
	dKLR(1,6,2,2,ttsf,gensf) =-hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2*CBA
	dKLR(6,1,1,1,ttsf,gensf) = dKLR(1,6,1,1,ttsf,gensf)
	dKLR(6,1,2,2,ttsf,gensf) = dKLR(1,6,2,2,ttsf,gensf)

	dKLR(2,6,1,1,ttsf,gensf) = EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*(-SA)/SB-Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*CA/CB)
	dKLR(2,6,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2*SBA
	dKLR(6,2,1,1,ttsf,gensf) = dKLR(2,6,1,1,ttsf,gensf)
	dKLR(6,2,2,2,ttsf,gensf) = dKLR(2,6,2,2,ttsf,gensf)

	dKLR(3,6,1,1,ttsf,gensf) =-cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MW2*TB)-Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2*TB) 
	dKLR(3,6,2,2,ttsf,gensf) = cI*hf_(ttsf,gensf)*hf_(ttsf+1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/sqrt2
	dKLR(6,3,1,1,ttsf,gensf) = dKLR(3,6,1,1,ttsf,gensf)
	dKLR(6,3,2,2,ttsf,gensf) = dKLR(3,6,2,2,ttsf,gensf)

	dKLR(4,6,1,1,ttsf,gensf) = cI*EL**2/(2d0*sqrt2*SW2)*(Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2-Mf2(ttsf+1,gensf)*(2*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))/MW2) 
	dKLR(4,6,2,2,ttsf,gensf) = (0d0,0d0)
	dKLR(6,4,1,1,ttsf,gensf) = dKLR(4,6,1,1,ttsf,gensf)
	dKLR(6,4,2,2,ttsf,gensf) = dKLR(4,6,2,2,ttsf,gensf)

	else

!h,h

	dKLR(1,1,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SA2                      
	dKLR(1,1,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SA2                 

!H,H
 
        dKLR(2,2,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CA2   
	dKLR(2,2,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CA2 

!P,P

        dKLR(3,3,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SB2    
	dKLR(3,3,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*SB2           

!G0,G0

        dKLR(4,4,1,1,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CB2 
	dKLR(4,4,2,2,ttsf,gensf) = -hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*CB2

!h,H

	dKLR(1,2,1,1,ttsf,gensf) = hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2A/2d0                   !LL-Kopplung
	dKLR(1,2,2,2,ttsf,gensf) = hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2A/2d0                   !RR-Kopplung
	dKLR(2,1,1,1,ttsf,gensf) = dKLR(1,2,1,1,ttsf,gensf)                                                              
	dKLR(2,1,2,2,ttsf,gensf) = dKLR(1,2,2,2,ttsf,gensf)                                                                          

!G0,P
      
	dKLR(3,4,1,1,ttsf,gensf) = hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2B/2d0 
	dKLR(3,4,2,2,ttsf,gensf) = hf_(ttsf,gensf)**2*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))*S2B/2d0 
	dKLR(4,3,1,1,ttsf,gensf) = dKLR(3,4,1,1,ttsf,gensf)
	dKLR(4,3,2,2,ttsf,gensf) = dKLR(3,4,2,2,ttsf,gensf)

!Hp-Hm

	dKLR(5,5,1,1,ttsf,gensf) = -2*EL**2/(MZ2*S2W**2)*Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/TB2
	dKLR(5,5,2,2,ttsf,gensf) = -2*EL**2*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MZ2*S2W**2)*TB2

!Gp-Gm

	dKLR(6,6,1,1,ttsf,gensf) = -2*EL**2/(MZ2*S2W**2)*Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))
	dKLR(6,6,2,2,ttsf,gensf) = -2*EL**2*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/(MZ2*S2W**2)

!Gp,Hm (Gm,Hp)

	dKLR(6,5,1,1,ttsf,gensf) =  - 2*EL**2/(MZ2*S2W**2)*Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/TB
	dKLR(6,5,2,2,ttsf,gensf) =  EL**2*TB/(2d0*SW2*MW2)*Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dKLR(5,6,1,1,ttsf,gensf) = dKLR(6,5,1,1,ttsf,gensf)
	dKLR(5,6,2,2,ttsf,gensf) = dKLR(6,5,2,2,ttsf,gensf)

!neutral Higgs,Hp (Hm)
!error in the HiggsHunter p.408, necessary replacements->see Couplings.F

	dKLR(1,5,1,1,ttsf,gensf) =-EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/MW2*CA*CB/SB2-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*(-SA)*SB/CB2) 
	dKLR(1,5,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2*SBA
	dKLR(5,1,1,1,ttsf,gensf) = dKLR(1,5,1,1,ttsf,gensf)
	dKLR(5,1,2,2,ttsf,gensf) = dKLR(1,5,2,2,ttsf,gensf)

	dKLR(2,5,1,1,ttsf,gensf) =-EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/MW2*SA*CB/SB2+Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*(-CA)*SB/CB2)
	dKLR(2,5,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2*CBA
	dKLR(5,2,1,1,ttsf,gensf) = dKLR(2,5,1,1,ttsf,gensf)
	dKLR(5,2,2,2,ttsf,gensf) = dKLR(2,5,2,2,ttsf,gensf)

	dKLR(3,5,1,1,ttsf,gensf) = cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/(MW2*TB2)+Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*TB2) 
	dKLR(3,5,2,2,ttsf,gensf) = (0d0,0d0)
	dKLR(5,3,1,1,ttsf,gensf) = dKLR(3,5,1,1,ttsf,gensf)
	dKLR(5,3,2,2,ttsf,gensf) = dKLR(3,5,2,2,ttsf,gensf)

	dKLR(4,5,1,1,ttsf,gensf) = cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/(MW2*TB)-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*TB) 
	dKLR(4,5,2,2,ttsf,gensf) = cI*hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2
	dKLR(5,4,1,1,ttsf,gensf) = dKLR(4,5,1,1,ttsf,gensf)
	dKLR(5,4,2,2,ttsf,gensf) = dKLR(4,5,2,2,ttsf,gensf)

!neutral Higgs,Gp (Gm)
!error in the HiggsHunter p.408, necessary replacements->see Couplings.F

	dKLR(1,6,1,1,ttsf,gensf) = EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/MW2*(-CA)/SB+Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*SA/CB) 
	dKLR(1,6,2,2,ttsf,gensf) =-hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2*CBA
	dKLR(6,1,1,1,ttsf,gensf) = dKLR(1,6,1,1,ttsf,gensf)
	dKLR(6,1,2,2,ttsf,gensf) = dKLR(1,6,2,2,ttsf,gensf)

	dKLR(2,6,1,1,ttsf,gensf) = EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/MW2*(-SA)/SB-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*CA/CB)
	dKLR(2,6,2,2,ttsf,gensf) = hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2*SBA
	dKLR(6,2,1,1,ttsf,gensf) = dKLR(2,6,1,1,ttsf,gensf)
	dKLR(6,2,2,2,ttsf,gensf) = dKLR(2,6,2,2,ttsf,gensf)

	dKLR(3,6,1,1,ttsf,gensf) = cI*EL**2/(2d0*sqrt2*SW2)*(-Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/(MW2*TB)-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2*TB) 
	dKLR(3,6,2,2,ttsf,gensf) =-cI*hf_(ttsf,gensf)*hf_(ttsf-1,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/sqrt2
	dKLR(6,3,1,1,ttsf,gensf) = dKLR(3,6,1,1,ttsf,gensf)
	dKLR(6,3,2,2,ttsf,gensf) = dKLR(3,6,2,2,ttsf,gensf)

	dKLR(4,6,1,1,ttsf,gensf) =-cI*EL**2/(2d0*sqrt2*SW2)*(Mf2(ttsf-1,gensf)*(2*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))/MW2-Mf2(ttsf,gensf)*(2*dMf(ttsf,gensf)/Mf(ttsf,gensf))/MW2) 
	dKLR(4,6,2,2,ttsf,gensf) = (0d0,0d0)
	dKLR(6,4,1,1,ttsf,gensf) = dKLR(4,6,1,1,ttsf,gensf)
	dKLR(6,4,2,2,ttsf,gensf) = dKLR(4,6,2,2,ttsf,gensf)

	endif
        
! ------------------------------ Vertex Squark - Squark - neutral Higgs - neutral Higgs ------------------------------

	do typex=1,4
          do typey = 1,4
	    do sfmeigenin=1,2
	       do sfmeigenout=1,2
	          do sfgaugeeigenin=1,2
	             do sfgaugeeigenout=1,2

    	 dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) +
     &   dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf,gensf)	               
   	            enddo

		dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) + 
     &   1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf,gensf) +               
     &                 dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf,gensf))
	            enddo
 	         enddo
	      enddo
       	   enddo
        enddo

! ------------------------------ Vertex Squark - Squark - charged Higgs - charged Higgs ------------------------------

	do typex=5,6
          do typey=5,6
	    do sfmeigenin=1,2
	       do sfmeigenout=1,2
	          do sfgaugeeigenin=1,2
	             do sfgaugeeigenout=1,2

    	 dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) +
     &   dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf,gensf)	               
   	            enddo

		dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf,gensf) + 
     &   1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf,gensf) +               
     &                 dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf,gensf))
	            enddo
 	         enddo
	      enddo
       	   enddo
        enddo

! ------------------------------ Vertex Squark - Squark - charged Higgs -neutral Higgs------------------------------
	do typex=5,6
          do typey=1,4
            do sfmeigenin=1,2
	      do sfmeigenout=1,2
! --- incoming up-type sfermion ---	
	      if ((ttsf.eq.1).or.(ttsf.eq.3)) then
	
	      do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
	            dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) +
     &              dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                   enddo
		    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf+1,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf+1,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf+1,gensf))
	         enddo

! --- incoming down-type sfermion ---
                    
                 else
                 do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
                    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) +
     &              dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                enddo
		    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf-1,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf-1,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf-1,gensf))
 		
	         enddo
                 endif

              enddo
            enddo
	  enddo
	enddo


! ------------------------------ Vertex Squark - Squark - neutral Higgs -charged Higgs------------------------------
	do typex=1,4
          do typey=5,6
            do sfmeigenin=1,2
	      do sfmeigenout=1,2
! --- incoming up-type sfermion ---	
	      if ((ttsf.eq.1).or.(ttsf.eq.3)) then
	
	      do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
	            dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) +
     &              dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                   enddo
		    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf+1,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf+1,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf+1,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf+1,gensf))
	         enddo

! --- incoming down-type sfermion ---
                    
                 else
                 do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
                    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) +
     &              dKLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     &              KLR(typex,typey,sfgaugeeigenin,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                enddo
		    dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) = dHHSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsf,ttsf-1,gensf) + 
     &              1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf-1,gensf))*HHSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsf,ttsf-1,gensf) +               
     &              dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*HHSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsf,ttsf-1,gensf))
 		
	         enddo
                 endif

              enddo
            enddo
	  enddo
	enddo	

	enddo
	enddo


!======================================================================================================================
!======================================  Vertex sFermion-sFermion-Vector-Vector =========================================
!======================================================================================================================
!see VVSfSf


        do typex=1,4                 !no flavor-mixing, theref_ore sfgen2 not needed
          do typey=1,4
            do sfmeigenin=1,2
              do sfmeigenout= 1,2
                do ttsfa=1,4
                  do gensf=1,3
                    do ttsfb=1,4

                   dVVSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = (0d0,0d0)

                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do ttsfa= 1, 4                    
          do ttsfb=1, 4
            do gensf=1, 3
              
        do sfmeigenin=1,2
          do sfmeigenout=1,2

!          if((ttsfa.eq.3).or.(ttsfa.eq.4)) then  

!g,g
       
!          dVVSfSf(1,1,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 0d0                                                                               !see Sparticles p. 89; EL=g2*SW, FeynArts MSSMQCD.dvi p.3              

!G,g

!	  dVVSfSf(1,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 0d0                                                                                !see  FeynArts MSSMQCD.dvi p.3 + NeuQ2qx_couplings 
!          dVVSfSf(2,1,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(1,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf)                                            

!g,Z

!	  dVVSfSf(1,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 2d0*EL*GS/(SW*CW)*(I3f_(ttsfa)*(dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf)+USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa,ttsfb)  !see FeynArts MSSMQCD.dvi p.3
!	  dVVSfSf(3,1,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(1,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf)   

!g,W+
!           if(ttsfa.eq.3) then
!            dVVSfSf(1,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = sqrt2*EL*GS/SW*(dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf)+USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf))*Krondelta(ttsfa+1,ttsfb)  !see FeynArts MSSMQCD.dvi p.4
!            dVVSfSf(4,1,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(1,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) 
!           else
!            dVVSfSf(1,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = sqrt2*EL*GS/SW*(dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf)+USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf))*Krondelta(ttsfa-1,ttsfb)  !see FeynArts MSSMQCD.dvi p.4
!            dVVSfSf(4,1,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(1,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) 
!           endif     
   
!          endif             

!GG

          dVVSfSf(2,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 0d0    !see Sparticles p. 512; EL=g2*SW, FeynArts MSSM.dvi p. 38

!G,Z
 
          dVVSfSf(2,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 2d0*EL**2/(SW*CW)*ef_(ttsfa)*I3f_(ttsfa)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf)+ USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa,ttsfb) !see FeynArts MSSM.dvi p. 38 
	  dVVSfSf(3,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) =  dVVSfSf(2,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf)      

!Z,Z
 
          dVVSfSf(3,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = 2d0*EL**2/(S2W**2)*(1d0 - 8*ef_(ttsfa)*I3f_(ttsfa)*SW2)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa,ttsfb) !see Feynarts MSSM.dvi p. 39 + ???Sparticles p. 513???      

!W+,W-

          dVVSfSf(4,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = EL**2/(2*SW2)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa,ttsfb) !see Feynarts MSSM.dvi p. 39 + ???Sparticles p. 513???       

!W+-,G/Z

          if((ttsfa.eq.1).or.(ttsfa.eq.3)) then

            dVVSfSf(2,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = EL**2/(3d0*sqrt2*SW)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa+1,ttsfb)  !see FeynArts MSSM.dvi p.39
            dVVSfSf(4,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(2,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) 


            dVVSfSf(3,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = -EL**2/(3d0*sqrt2*CW)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa+1,ttsfb)  !see FeynArts MSSM.dvi p.39
            dVVSfSf(4,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(3,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf)  

          else

            dVVSfSf(2,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = EL**2/(3d0*sqrt2*SW)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa-1,ttsfb)  !see FeynArts MSSM.dvi p.39
            dVVSfSf(4,2,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(2,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) 


            dVVSfSf(3,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = -EL**2/(3d0*sqrt2*CW)*( dUSf(sfmeigenin,1,ttsfa,gensf)*USf(sfmeigenout,1,ttsfb,gensf) + USf(sfmeigenin,1,ttsfa,gensf)*dUSf(sfmeigenout,1,ttsfb,gensf) )*Krondelta(ttsfa-1,ttsfb)  !see FeynArts MSSM.dvi p.39
            dVVSfSf(4,3,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(3,4,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf)  

          endif 


          do typex = 2,4                      !no Gluon-Counterterm needed yet => Sum from 2-4 enough so far
            do typey = 2,4
              do sfgaugeeigenin = 1,2

          dVVSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) = dVVSfSf(typex,typey,sfmeigenin,sfmeigenout,ttsfa,ttsfb,gensf) + 
     &                                                                    1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsfb,gensf))*VVSfSf(typex,typey,sfmeigenin,sfgaugeeigenin,ttsfa,ttsfb,gensf) +               
     &                                                                    dZSf(sfgaugeeigenin,sfmeigenin,ttsfa,gensf)*VVSfSf(typex,typey,sfgaugeeigenin,sfmeigenout,ttsfa,ttsfb,gensf))
             
              enddo
            enddo
          enddo
               
                enddo
              enddo
            enddo
          enddo
        enddo

 	end
        
! =======================================
! ======= Propagator counterterms ======= 
! =======================================

! Squarks
	subroutine QQ2xx_CalcPropSquark(myP2,type,gen)

#include "QQ2xx_Model.h"
#include "QQ2xx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer type,gen,i,j
	double precision myP2, Krondelta

	do i = 1, 2
	  do j = 1, 2

	    Prop(i,j,type,gen) = 1/2d0*dZSf(i,j,type,gen)*(myP2 - MSf2(i,type,gen)) + 1/2d0*dZSf(j,i,type,gen)*(myP2 - MSf2(j,type,gen)) - Krondelta(i,j)*dMSf2(i,type,gen)
      

!            print*, "Prop((i,j,type,gen)", Prop(i,j,type,gen),type, gen

	  enddo
 	enddo

	end


! Quarks
	subroutine QQ2xx_CalcPropQuark(type,gen)

#include "QQ2xx_Model.h"
#include "QQ2xx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer type,gen

! PropVR/L according to PL/PR before kslash...
! dZfl, ... etc. follows convention like ChiChi2QQ code...

	    PropVR(type,gen) = 1/2d0*(dZfL(type,gen) + dconjg(dZfL(type,gen)))
	    PropVL(type,gen) = 1/2d0*(dZfR(type,gen) + dconjg(dZfR(type,gen)))
	    PropSL(type,gen) = (-1)*1/2d0*(Mf(type,gen)*dZfL(type,gen) + Mf(type,gen)*dconjg(dZfR(type,gen)))- dMf(type,gen)
	    PropSR(type,gen) = (-1)*1/2d0*(Mf(type,gen)*dZfR(type,gen) + Mf(type,gen)*dconjg(dZfL(type,gen)))- dMf(type,gen)

c             print*, "PropVR(type,gen)", PropVR(type,gen)
c             print*, "PropVL(type,gen)", PropVL(type,gen)
c             print*, "PropSR(type,gen)", PropSR(type,gen)
c             print*, "PropSL(type,gen)", PropSL(type,gen)

	end
