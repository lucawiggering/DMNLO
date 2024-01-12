! =====================================================================
	subroutine stst2QQ_SetCounterTerms()
	implicit none

#include "stst2QQ_Model.h"
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_GenCouplings.h"
#include "stst2QQ_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer a,b,i,j,k,l,m,t,g,n,isfin,isfout,ittin,ittout,igenin,igenout,o,p,tt,gen
	integer sfmeigenout,eigensf,ttsf,gensf, igensf,ittsf
	integer sfout, sftypeout, ftypein,sfgenout, fgenin
	integer xtype, typex, ttsfb, ttsfa,sfgaugeeigenout,sfgaugeeigenin,gensfa,gensfb,ytype
	integer eigenchi
	double precision dCLnew, CL(2,2,4,3,4), dCL(2,2,4,3,4),dBRtest,dBLtest
	double precision mbmass,mtmass, mtit,CLR11, CLR22, CL1, CL2
        double complex DLred(6),DRred(6), doof

! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf

!                 For Debugging only:
!                 DeltaGs = 0d0
!                 dZG = 0d0
!                 dZGlL= 0d0
!                 dZGlR= 0d0
! !                 dZfL=0d0
! !                 dZfR=0d0
! ! 	   dZGlR = 0d0
! ! 	   dZGlL = 0d0
!             dMGl = 0d0

!************************
! Coupling counterterms *
!************************


! 
! ******************************************************** 
! c Vertex Gluino  - Squark (outgoing) - Quark (incoming) *
! ********************************************************
! c e.g. GluL/R(sfmeigenin,sftypein,ftypeout,sfgenin,fgenout), couplings multiplied by -i
! 
           do sfout=1, 2                                                                              
             do sftypeout= 1, 4
               do ftypein = 1, 4
                 do sfgenout= 1, 3
                    do fgenin= 1, 3

                    if((sftypeout.gt.2).and.(ftypein.gt.2)) then

			 dGluL(sfout,sftypeout,ftypein,sfgenout,fgenin) = - GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,1,sftypeout,sfgenout)+dUSf(sfout,1,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluL(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluL(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfL(ftypein,fgenin)*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin)+dZGlL*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
     
                        dGluR(sfout,sftypeout,ftypein,sfgenout,fgenin) = GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,2,sftypeout,sfgenout)+dUSf(sfout,2,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluR(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluR(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfR(ftypein,fgenin)*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin)+dZGlR*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
     
			dGluLt1(sfout,sftypeout,ftypein,sfgenout,fgenin) = GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,2,sftypeout,sfgenout)+dUSf(sfout,2,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluR(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluR(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfR(ftypein,fgenin)*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin)+dZGlL*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
			dGluRt1(sfout,sftypeout,ftypein,sfgenout,fgenin) = - GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,1,sftypeout,sfgenout)+dUSf(sfout,1,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluL(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluL(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfL(ftypein,fgenin)*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin)+dZGlR*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
			
			dGluRt2(sfout,sftypeout,ftypein,sfgenout,fgenin) = GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,2,sftypeout,sfgenout)+dUSf(sfout,2,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluR(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluR(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfR(ftypein,fgenin)*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin)+dconjg(dZGlL)*GluR(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
			dGluLt2(sfout,sftypeout,ftypein,sfgenout,fgenin) =  - GS*sqrt2*Krondelta(sftypeout,ftypein)*Krondelta(sfgenout,fgenin)*(DeltaGs*USf(sfout,1,sftypeout,sfgenout)+dUSf(sfout,1,sftypeout,sfgenout))+
     & 							1d0/2d0*(dZSf(1,sfout,sftypeout,sfgenout)*GluL(1,sftypeout,ftypein,sfgenout,fgenin)+dZSf(2,sfout,sftypeout,sfgenout)*GluL(2,sftypeout,ftypein,sfgenout,fgenin)+
     &							dZfL(ftypein,fgenin)*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin)+dconjg(dZGlR)*GluL(sfout,sftypeout,ftypein,sfgenout,fgenin))
     
                     else

                        dGluL(sfout,sftypeout,ftypein,sfgenout,fgenin) = (0d0,0d0)
                        dGluR(sfout,sftypeout,ftypein,sfgenout,fgenin) = (0d0,0d0)
                        
                        dGluRt1(sfout,sftypeout,ftypein,sfgenout,fgenin) = dGluL(sfout,sftypeout,ftypein,sfgenout,fgenin)
                        dGluLt1(sfout,sftypeout,ftypein,sfgenout,fgenin) = dGluR(sfout,sftypeout,ftypein,sfgenout,fgenin)

                        dGluRt2(sfout,sftypeout,ftypein,sfgenout,fgenin) = dGluR(sfout,sftypeout,ftypein,sfgenout,fgenin)
                        dGluLt2(sfout,sftypeout,ftypein,sfgenout,fgenin) = dGluL(sfout,sftypeout,ftypein,sfgenout,fgenin)
                        
                     endif     

                    enddo
                 enddo
              enddo
            enddo
          enddo
	  
! ===== Neu-Sf-f -- A-coupling =====
! 
! !  e.g. AL/R(chimeigen, sfmeigen, sftype, sfgen)
! 	 
         do eigensf = 1, 2
 	   do ttsf = 1, 4
	    do k=1,4
 	      do gensf = 1, 3
		do m=1,3
!  	      AL/R(chimeigen,sfmeigenin,sftypein,ftypeout,sfgenin,fgenout)

! ------------------------------ Vertex Neutralino - Sfermion (in) - Fermion (out) ------------------------------

 	         do eigenchi = 1, 4

! --- incoming up-type sfermion ---

 	               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
 	               
 	               	  dAR(eigenchi,eigensf,ttsf,k,gensf,m) = (- EL/(SW*sqrt2*MW*SB)*ZNeu(eigenchi,4)*(dMf(k,m)*USf(eigensf,2,ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)*Mf(k,m)) -
     &                    sqrt2*EL/SW*((ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1) +  I3f(ttsf)*ZNeu(eigenchi,2))*(dUSf(eigensf,1,ttsf,gensf)))*Krondelta(ttsf,k)*Krondelta(gensf,m)+
     &                    AR(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AR(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                    AR(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfL(ttsf,gensf))
     
                          dAL(eigenchi,eigensf,ttsf,k,gensf,m) =  (- EL/(SW*sqrt2*MW*SB)*ZNeuC(eigenchi,4)*(dMf(k,m)*USf(eigensf,1,ttsf,gensf) + dUSf(eigensf,1,ttsf,gensf)*Mf(k,m)) + 
     &                    sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*(dUSf(eigensf,2,ttsf,gensf)))*Krondelta(ttsf,k)*Krondelta(gensf,m)+
     & 	                  AL(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AL(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                    AL(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfR(ttsf,gensf))	
     
     
			  dALt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAL(eigenchi,eigensf,ttsf,k,gensf,m)
			  dARt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAR(eigenchi,eigensf,ttsf,k,gensf,m)
			  
			  dALt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAR(eigenchi,eigensf,ttsf,k,gensf,m))
			  dARt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAL(eigenchi,eigensf,ttsf,k,gensf,m))
  
    
! ! --- incoming down-type sfermion---     

      		       else
!       		       
                           dAR(eigenchi,eigensf,ttsf,k,gensf,m) =  (- EL/(SW*sqrt2*MW*CB)*ZNeu(eigenchi,3)*(dMf(k,m)*USf(eigensf,2,ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)*Mf(k,m)) -
     &                     sqrt2*EL/SW*(I3f(ttsf)*ZNeu(eigenchi,2) + (ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1))*(dUSf(eigensf,1,ttsf,gensf)))*Krondelta(ttsf,k)*Krondelta(gensf,m)+
     & 	                   AR(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AR(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                     AR(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfL(ttsf,gensf))
! 
                           dAL(eigenchi,eigensf,ttsf,k,gensf,m) =  (- EL/(SW*sqrt2*MW*CB)*ZNeuC(eigenchi,3)*(dMf(k,m) *USf(eigensf,1,ttsf,gensf) + dUSf(eigensf,1,ttsf,gensf)*Mf(k,m)) + 
     &                     sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*(dUSf(eigensf,2,ttsf,gensf)))*Krondelta(ttsf,k)*Krondelta(gensf,m)+
     & 	                   AL(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AL(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                     AL(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfR(ttsf,gensf))
     
     
			  dALt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAL(eigenchi,eigensf,ttsf,k,gensf,m)
			  dARt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAR(eigenchi,eigensf,ttsf,k,gensf,m)
			  
			  dALt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAR(eigenchi,eigensf,ttsf,k,gensf,m))
			  dARt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAL(eigenchi,eigensf,ttsf,k,gensf,m))
! 			  			  

	               endif
	          enddo

	          

! ------------------------------ Vertex Chargino - Sfermion (in) - Fermion (out) ------------------------------

	          do eigenchi= 5, 6   

! --- incoming up-type sfermion ---	          

	               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
            
 	         	dAL(eigenchi,eigensf,ttsf,k,gensf,m) = EL*Mf(ttsf+1,gensf)/(SW*sqrt2*MW*CB)*dconjg(UCha(eigenchi-4,2))*
     &          USf(eigensf,1,ttsf,gensf)*Krondelta(ttsf,k-1)*Krondelta(gensf,m)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf)+
     &          dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + AL(eigenchi,1,ttsf,k,gensf,m)*1d0/2d0*dZSf(1,eigensf,ttsf,gensf)+ 
     &          AL(eigenchi,2,ttsf,k,gensf,m)*1d0/2d0*dZSf(2,eigensf,ttsf,gensf) + AL(eigenchi,eigensf,ttsf,k,gensf,m)*1d0/2d0*dconjg(dZfR(k,gensf))
 	         	
 	         	dAR(eigenchi,eigensf,ttsf,k,gensf,m) = (-EL/SW*VCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf)*
     &          (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)*
     &          USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf)+dUSf(eigensf,2,ttsf,gensf)/
     &          USf(eigensf,2,ttsf,gensf)))*Krondelta(ttsf,k-1)*Krondelta(gensf,m) + AR(eigenchi,1,ttsf,k,gensf,m)*
     &          1d0/2d0*dZSf(1,eigensf,ttsf,gensf)+AR(eigenchi,2,ttsf,k,gensf,m)*1d0/2d0*dZSf(2,eigensf,ttsf,gensf)+
     &          AR(eigenchi,eigensf,ttsf,k,gensf,m)*1d0/2d0*dconjg(dZfL(k,gensf))
!  	         	
! 
! 
! ! --- incoming down-type sfermion --- 	 
! 
 	         	else


!  	         	dAL(eigenchi,eigensf,ttsf,k,gensf,m) = EL*Mf(ttsf-1,gensf)/(SW*sqrt2*MW*SB)*VChaC(eigenchi-4,2)*USf(eigensf,1,ttsf,gensf)*Krondelta(ttsf,k+1)*Krondelta(gensf,m)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))) + 
!      &                  AL(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
!      &                  AL(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
!      &                  AL(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfR(ttsf-1,gensf))

 	         	dAL(eigenchi,eigensf,ttsf,k,gensf,m) = EL*Mf(ttsf-1,gensf)/(SW*sqrt2*MW*SB)*dconjg(VCha(eigenchi-4,2))*
     &          USf(eigensf,1,ttsf,gensf)*Krondelta(ttsf,k+1)*Krondelta(gensf,m)*
     &          (dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf)+dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))+ 
     &                  AL(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfR(ttsf-1,gensf))

 	         	dAR(eigenchi,eigensf,ttsf,k,gensf,m) = -EL/SW*UCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf)*
     &          (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)*
     &          USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf)) + 
     &                  AR(eigenchi,1,ttsf,k,gensf,m)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,2,ttsf,k,gensf,m)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,eigensf,ttsf,k,gensf,m)*1/2d0*dconjg(dZfL(ttsf-1,gensf))
     
 	         	endif
 	         	
              dALt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAL(eigenchi,eigensf,ttsf,k,gensf,m)
			  dARt1(eigenchi,eigensf,ttsf,k,gensf,m) = dAR(eigenchi,eigensf,ttsf,k,gensf,m)
			  
			  dALt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAR(eigenchi,eigensf,ttsf,k,gensf,m))
			  dARt2(eigenchi,eigensf,ttsf,k,gensf,m) = dconjg(dAL(eigenchi,eigensf,ttsf,k,gensf,m))
			  
	          enddo
 	      enddo
	   enddo
	enddo
	enddo
	enddo


	end
