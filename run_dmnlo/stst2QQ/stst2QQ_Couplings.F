	subroutine stst2qq_SetCouplings
	implicit none

#include "stst2QQ_Model.h"
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_GenCouplings.h"
#include "stst2QQ_GenCounterterms.h"

	integer i, j, k, l, m, n, xtype, ytype, a, b, o, p, tt1, tt2, gen1, gen2,sfmeigenout1,sfmeigenin1

! functions from DMNLO_Kinematics.F
      double precision Krondelta, ef, I3f, hf
      double precision DmbResults(5)
      external DMNLO_CalcDeltaMb

! ********************************
!  Fermion-Vector boson couplings*
! ********************************
!e.g. BL(i,f3type,f4type,f3gen,f4gen), no generation-mixing, couplings multiplied by -i

!1=Gluon, 2= Photon, 3= Z-Boson, 4= W-Boson

       do i = 1, 4 
         do j = 1, 4
            do m = 1, 3
               do n = 1, 3
         
 
                if((i.gt.2).and.(j.gt.2)) then

                  BL(1,i,j,m,n) = -GS*Krondelta(i,j)*Krondelta(m,n)
                  BR(1,i,j,m,n) = -GS*Krondelta(i,j)*Krondelta(m,n)
 
                else

                  BL(1,i,j,m,n) = (0d0,0d0)
                  BR(1,i,j,m,n) = (0d0,0d0)

                endif
                                
                BL(2,i,j,m,n) = -ef(i)*EL*Krondelta(i,j)*Krondelta(m,n)
                BR(2,i,j,m,n) = -ef(i)*EL*Krondelta(i,j)*Krondelta(m,n)                

                BV(i) = I3f(i) - 2d0*SW**2*ef(i)
                BA(i) = I3f(i)

                BL(3,i,j,m,n) = - (EL/(2d0*CW*SW))*(BV(i) + BA(i))*Krondelta(i,j)*Krondelta(m,n)
                BR(3,i,j,m,n) = - (EL/(2d0*CW*SW))*(BV(i) - BA(i))*Krondelta(i,j)*Krondelta(m,n)


                if((i.eq.1).or.(i.eq.3)) then               !incoming up type (->W+ -radiation)

                BL(4,i,j,m,n) = - EL/(sqrt2*SW)*Krondelta(i,j-1)*Krondelta(m,n)

                else                                        !incoming down type (->W- -radiation)

                BL(4,i,j,m,n) = - EL/(sqrt2*SW)*Krondelta(i,j+1)*Krondelta(m,n)

                endif

                BR(4,i,j,m,n) = (0d0,0d0)

              enddo
            enddo
          enddo
	enddo
	
	                
!         write(*,*)'BL(3,3,3,3)',BL(3,3,3,3,3)
!         write(*,*)'BR(3,3,3,3)',BR(3,3,3,3,3)

!**************************	
! Fermion-Higgs couplings *
!**************************      
!e.g. DL/R(i,f3type,f4type,f3gen,f4gen), couplings multiplied by -i

! =====================================================================================================	
! ==================================== Vertex Quark - Quark - Higgs ===================================
! =====================================================================================================

! e.g.  FL/R(xtype,pstype,psgen) 
! 	      

!  
! 	write(*,*)'Kop check Jul',-EL*Mf(3,3)*CA/(SW*2d0*MW*SB) 
! 	write(*,*)'Kop check me',-hf(3,3)*CA/sqrt2
        do i = 1, 4
           do j = 1, 4
              do m = 1, 3
                 do n = 1, 3
!                  DL(1,i,j,m,n) = 1d0
!                  DR(1,i,j,m,n) = 1d0
!            
!--- incoming up-type fermion ---           
!  hf = EL*Mf(type,gen)/(sqrt2*MW*SW*SB) 
	
           
              if ((i.eq.1).or.(i.eq.3)) then
! h              
                  DL(1,i,j,m,n) = -hf(i,m)*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)                
                  DR(1,i,j,m,n) =  DL(1,i,j,m,n)
! H                  
                  DL(2,i,j,m,n) = -hf(i,m)*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(2,i,j,m,n) =  DL(2,i,j,m,n)
! H3                  
                  DL(3,i,j,m,n) = -cI*hf(i,m)*CB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
                  DR(3,i,j,m,n) = -DL(3,i,j,m,n)
! G0                  
                  DL(4,i,j,m,n) = -cI*hf(i,m)*SB/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(4,i,j,m,n) = -DL(4,i,j,m,n)
! H+(out)                
		  DL(5,i,j,m,n) =  hf(i+1,m)*SB*Krondelta(i,j-1)*Krondelta(m,n)                  
		  DR(5,i,j,m,n) =  hf(i,m)*CB*Krondelta(i,j-1)*Krondelta(m,n)
! G+(out)                 
	          DL(6,i,j,m,n) = -hf(i+1,m)*CB*Krondelta(i,j-1)*Krondelta(m,n)
		  DR(6,i,j,m,n) = hf(i,m)*SB*Krondelta(i,j-1)*Krondelta(m,n)
 
!                   
! --- incoming down-type fermion --- 

 	             else
! h 	
 	          DL(1,i,j,m,n) = hf(i,m)*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)   
 	          DR(1,i,j,m,n) = DL(1,i,j,m,n)
! H 	          
 	          DL(2,i,j,m,n) = -hf(i,m)*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
 	          DR(2,i,j,m,n) =  DL(2,i,j,m,n)
! H3 	          
 	          DL(3,i,j,m,n) = -cI*hf(i,m)*SB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
 	          DR(3,i,j,m,n) = -DL(3,i,j,m,n)
! G0 	          
 	          DL(4,i,j,m,n) = cI*hf(i,m)*CB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
 	          DR(4,i,j,m,n) = -DL(4,i,j,m,n)
! H-(out) 	          
 	          DL(5,i,j,m,n) = hf(i-1,m)*CB*Krondelta(i,j+1)*Krondelta(m,n)
		  DR(5,i,j,m,n) = hf(i,m)*SB*Krondelta(i,j+1)*Krondelta(m,n) 
! G-(out) 	          
		  DL(6,i,j,m,n) = hf(i-1,m)*SB*Krondelta(i,j+1)*Krondelta(m,n)                  
		  DR(6,i,j,m,n) = -hf(i,m)*CB*Krondelta(i,j+1)*Krondelta(m,n)
 
 	           	          
 	          endif

                enddo
             enddo
 	   enddo
 	enddo
 	
 
	 
!
!        
       
!******************************************************** 
! Vertex Gluino  - Squark (outgoing) - Quark (incoming) *
!********************************************************
!  e.g. GluL/R(sfmeigenout,sftypeout,ftypein,sfgenout,fgenin), couplings multiplied by -i

           do a=1, 2                                                                              !taken from NeuQ2qx + Sparticles p. 88 
             do b= 1, 4
               do m= 1, 4
                 do i= 1, 3
                    do n= 1, 3

                     if((b.gt.2).and.(m.gt.2)) then

                        GluL(a,b,m,i,n) = - GS*sqrt2*USf(a,1,b,i)*Krondelta(b,m)*Krondelta(i,n)
                        GluR(a,b,m,i,n) =   GS*sqrt2*USf(a,2,b,i)*Krondelta(b,m)*Krondelta(i,n)
                        ! GluL/R(sfmeigenout,sftypeout,ftypein,sfgenout,fgenin)
                        GluLt2(a,b,m,i,n) = GluL(a,b,m,i,n)
                        GluRt2(a,b,m,i,n) = GluR(a,b,m,i,n)
                        ! GluL/R(sfmeigenin,sftypein,ftypeou,sfgenin,fgenout)
                        GluLt1(a,b,m,i,n) = GluR(a,b,m,i,n)
                        GluRt1(a,b,m,i,n) = GluL(a,b,m,i,n)

                     else

                        GluL(a,b,m,i,n) = (0d0,0d0)
                        GluR(a,b,m,i,n) = (0d0,0d0)

                     endif     


                    enddo
                 enddo
              enddo
            enddo
          enddo
   
   
! ***********************************************
!Sfermion - Sfermion - Gluon / Photon-couplings *             
!************************************************
! e.g. C(xtype,sftype1) 
! e.g. VSfSf(xtype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2) / general structure (s-channel):VSfSf*(p1-p2), couplings multiplied by -i

!1=Gluon (=g), 2=Gamma (=G)

	C(1,1) =  (0d0,0d0)
	C(1,2) =  (0d0,0d0)


	C(1,3) = -GS
	C(1,4) = -GS

	C(2,1) = -ef(1)*EL
	C(2,2) = -ef(2)*EL
	C(2,3) = -ef(3)*EL
	C(2,4) = -ef(4)*EL

        do i = 1, 2
          do j = 1, 2   
            do k = 1, 4
              do l = 1, 4
                do m = 1, 3
                   do n = 1, 3

                       VSfSf(1,i,j,k,l,m,n) = C(1,k)*Krondelta(i,j)*Krondelta(k,l)*Krondelta(m,n)                     
                       VSfSf(2,i,j,k,l,m,n) = C(2,k)*Krondelta(i,j)*Krondelta(k,l)*Krondelta(m,n)                      

                       VSfSft(1,i,j,k,l,m,n) = VSfSf(1,i,j,k,l,m,n)                                                    
                       VSfSfu(1,i,j,k,l,m,n) = VSfSf(1,i,j,k,l,m,n)                                                    
                       VSfSft(2,i,j,k,l,m,n) = VSfSf(2,i,j,k,l,m,n)                                                   
                       VSfSfu(2,i,j,k,l,m,n) = VSfSf(2,i,j,k,l,m,n)                                                   

                  enddo
                enddo
              enddo
       	    enddo
       	  enddo
      	enddo
  
        



! *****************************************  
! Sfermion - Sfermion - W / Z - couplings *              
! *****************************************
! e.g. VSfSf(xtype,sfmeigenout,sfmeigenin,sftypeout,sftypein,sfgenout,sfgenin) / general structure (s-channel):VSfSf*(p1-p2), couplings multiplied by -i

!3=Z, 4=W+ (W-)

        do i = 1, 2  
          do j = 1, 2    
            do k = 1, 4     
              do l = 1, 4
                do m = 1, 3
                  do n = 1, 3                                   
     
     
           VSfSf(3,i,j,k,l,m,n)=-EL/(SW*CW)*(I3f(k)*USf(j,1,l,n)*USf(i,1,k,m)-SW2*ef(k)*Krondelta(j,i))*Krondelta(k,l)*Krondelta(m,n)  !Z-Boson
           
          
           
           VSfSft(3,i,j,k,l,m,n) = VSfSf(3,i,j,k,l,m,n) 
           VSfSfu(3,i,j,k,l,m,n) = VSfSf(3,i,j,k,l,m,n)        
           
     
! --- incoming up-type sfermion ---
               
           if ((l.eq.1).or.(l.eq.3)) then                                                                                                !W+/- -Boson

    	    VSfSf(4,i,j,k,l,m,n) = - EL/(sqrt2*SW)*USf(j,1,l,n)*USf(i,1,k,m)*Krondelta(l+1,k)*Krondelta(m,n)                          


! --- incoming down-type sfermion ---
               
           else
    	       
    	    VSfSf(4,i,j,k,l,m,n) = - EL/(sqrt2*SW)*USf(j,1,l,n)*USf(i,1,k,m)*Krondelta(l-1,k)*Krondelta(m,n)         

                                
           endif
               

                    enddo
                  enddo 
                enddo
 	      enddo
 	    enddo
	  enddo


!*************************************
! Gaugino-Fermion-Sfermion couplings *
!*************************************
! e.g. AL/R(chimeigen,sfmeigenin,sftypein,ftypeout,sfgenin,fgenout), couplings multiplied by -i
	 
 	 
        do a = 1, 2
 	   do b = 1, 4
             do m= 1, 4
 	       do i = 1, 3
 	         do n =1, 3
 	      
! ------------------------------ Vertex Neutralino - Sfermion (in) - Fermion (out) ------------------------------

 	         do j = 1, 4  !taken from Sparticles p. 517 with GL->AR,GR->AL etc. (following the sfermion(in)-diagramm + ZneuC can be imaginary=> keep track of hermitian conjugation (see Sparticles p. 517))

! --- incoming up-type sfermion ---

 	               if ((b.eq.1).or.(b.eq.3)) then
 	               
 	                  AR(j,a,b,m,i,n) = (- EL*Mf(b,i)/(SW*sqrt2*MW*SB)*ZNeu(j,4)*USf(a,2,b,i) -
     &                    sqrt2*EL/SW*((ef(b) - I3f(b))*SW/CW*ZNeu(j,1) +  I3f(b)*ZNeu(j,2))*USf(a,1,b,i))*Krondelta(b,m)*Krondelta(i,n)
     
                      AL(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*SB)*ZNeuC(j,4)*USf(a,1,b,i) + 
     &                    sqrt2*EL*ef(b)/CW*ZNeuC(j,1)*USf(a,2,b,i))*Krondelta(b,m)*Krondelta(i,n)	
! 			 
			  
	  
			  ARt1(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ARt2(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ALt2(j,a,b,m,i,n) = AR(j,a,b,m,i,n)


! --- incoming down-type sfermion---     

      		       else
      		       
                          AR(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*CB)*ZNeu(j,3)*USf(a,2,b,i) -
     &                    sqrt2*EL/SW*(I3f(b)*ZNeu(j,2) + (ef(b) - I3f(b))*SW/CW*ZNeu(j,1))*USf(a,1,b,i))*Krondelta(b,m)*Krondelta(i,n)

                          AL(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*CB)*ZNeuC(j,3)*USf(a,1,b,i) + 
     &                    sqrt2*EL*ef(b)/CW*ZNeuC(j,1)*USf(a,2,b,i))*Krondelta(b,m)*Krondelta(i,n)
     
     
			  ARt1(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ARt2(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ALt2(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
                                 
	               endif
	          enddo

	          
! ------------------------------ Vertex Chargino - Sfermion (in) - Fermion (out) ------------------------------                                                         !differs from couplings implemented in CalcHEP->not comparable

	          do j= 5, 6   

          

	          if ((b.eq.1).or.(b.eq.3)) then
                                   

! --- incoming up-type sfermion ---	
	         	  AL(j,a,b,m,i,n) =  EL*Mf(b+1,i)/(SW*sqrt2*MW*CB)*dconjg(UCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m-1)*Krondelta(i,n)                             != E_R in Sparticles p. 515 , no t_i-e-coupl. etc. and no generation-mixing                                                                                                                                                                        
 	         	  AR(j,a,b,m,i,n) = (-EL/SW*VCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)/(SW*sqrt2*MW*SB)*VCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m-1)*Krondelta(i,n)    != D_L in Sparticles p. 515
! --- incoming up-type sfermion ---	     
			  ARt1(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
! --- outgoing up-type sfermion ---	
			  ARt2(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ALt2(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  

                        
 	         	else


! --- incoming down-type sfermion --- 
 	         	  AL(j,a,b,m,i,n) = EL*Mf(b-1,i)/(SW*sqrt2*MW*SB)*dconjg(VCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m+1)*Krondelta(i,n)                              != F_R in Sparticles p. 515
 	         	  AR(j,a,b,m,i,n) = (-EL/SW*UCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)/(SW*sqrt2*MW*CB)*UCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m+1)*Krondelta(i,n)    != C_L in Sparticles p. 515 
 ! --- incoming down-type sfermion ---      
			  ARt1(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
! --- outgoing down-type sfermion ---
			  ARt2(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  ALt2(j,a,b,m,i,n) = AR(j,a,b,m,i,n)	         	
 	         	endif
	          enddo

                  enddo
                enddo
 	      enddo
	   enddo
        enddo
        
 
 
 
!  
! ** Delta_Mb & the resummation
	  call DMNLO_CalcDeltaMb(DmbResults,3)
	
	! Default Setting is zero
        do a = 1, 2
 	     do b = 1, 4
          do m= 1, 4
 	       do i = 1, 3
 	        do n =1, 3
 ! --- incoming down-type sfermion ---      
			ARt1ho(j,a,b,m,i,n) = 0d0
			ALt1ho(j,a,b,m,i,n) = 0d0
            enddo
           enddo
          enddo
         enddo
        enddo
        
        do a = 1, 2
 	     do b = 1, 4
          do m= 1, 4
 	       do i = 1, 3
 	        do n =1, 3
  
 	         
 	          do j = 1, 4  !taken from Sparticles p. 517 with GL->AR,GR->AL etc. (following the sfermion(in)-diagramm + ZneuC can be imaginary=> keep track of hermitian conjugation (see Sparticles p. 517))

! --- incoming up-type sfermion ---[

 	               if ((b.eq.1).or.(b.eq.3)) then
 
			  ARt1ho(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1ho(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			  
! --- incoming down-type sfermion---     

      		       else
      		       
			  ARt1ho(j,a,b,m,i,n) = AR(j,a,b,m,i,n)
			  ALt1ho(j,a,b,m,i,n) = AL(j,a,b,m,i,n)
			
	               endif
	          enddo
	          
            do j= 5, 6 
! 
 	        if(b.eq.3) then
 ! --- incoming down-type sfermion ---      
			ARt1ho(j,a,b,m,i,n) = (-EL/SW*VCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)/(SW*sqrt2*MW*SB)*VCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m-1)*Krondelta(i,n)
			ALt1ho(j,a,b,m,i,n) = EL*Mf(b+1,i)*(1d0/(1d0 + DeltaMb) + DeltaMbQCD)/(SW*sqrt2*MW*CB)*dconjg(UCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m-1)*Krondelta(i,n)  
			elseif(b.eq.4) then
! --- outgoing down-type sfermion ---
			ARt1ho(j,a,b,m,i,n) = (-EL/SW*UCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)*(1d0/(1d0 + DeltaMb) + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m+1)*Krondelta(i,n) 
			ALt1ho(j,a,b,m,i,n) =  EL*Mf(b-1,i)/(SW*sqrt2*MW*SB)*dconjg(VCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m+1)*Krondelta(i,n)
			endif
			
			
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        
        
        
    
          
!********************************************
! sFermion-sFermion-Gluon-Gluon-couplings *
!********************************************

	do i=1,2
	  do j=1,2
	    do tt1=1,4
	      do tt2=1,4
	        do gen1=1,3
		  do gen2=1,3
		SfSfgg(i,tt1,gen1,j,tt2,gen2) = 0.0
		  enddo
		enddo
	      enddo
	    enddo
	  enddo
	enddo

	do i=1,2
	  do j=1,2
	    do tt1=1,4
	      do tt2=1,4
	        do gen1=1,3
		  do gen2=1,3
		SfSfgg(i,tt1,gen1,j,tt2,gen2) = GS*GS*Krondelta(i,j)*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)
		  enddo
		enddo
	      enddo
	    enddo
	  enddo
	enddo          

	
	         
!********************************************
! sFermion-sFermion-Gluon-Vectorboson-couplings *
!********************************************
! 
	do k=1,4
	do i=1,2
	  do j=1,2
	    do tt1=1,4
	      do tt2=1,4
	        do gen1=1,3
		  do gen2=1,3
		SfSfVV(k,i,tt1,gen1,j,tt2,gen2) = 0.0
		  enddo
		enddo
	      enddo
	    enddo
	  enddo
	enddo
	enddo

	do i=1,2
	  do j=1,2
	    do tt1=1,4
	      do tt2=1,4
	        do gen1=1,3
		  do gen2=1,3
		  SfSfVV(1,i,tt1,gen1,j,tt2,gen2) = SfSfgg(i,tt1,gen1,j,tt2,gen2) 
		  SfSfVV(2,i,tt1,gen1,j,tt2,gen2) = 2d0*EL*GS*ef(itt1)*Krondelta(i,j)*Krondelta(itt1,itt2)*Krondelta(igen1,igen2)    
		  SfSfVV(3,i,tt1,gen1,j,tt2,gen2) = 2d0*EL*GS/(SW*CW)*(I3f(itt1)*USf(i,1,itt1,igen1)*USf(j,1,itt2,igen2)-ef(itt1)*SW2*Krondelta(i,j))*Krondelta(itt1,itt2)*Krondelta(igen1,igen2)
		  enddo
		enddo
	      enddo
	    enddo
	  enddo
	enddo          

! ===================================================================================================== 
! =============================== Vertex Squark - Squark - Squark - Squark ============================
! =====================================================================================================
	

	do sfmeigenin1=1,2
	   do sfmeigenout1=1,2
	      do tt1 = 1, 4
	      do tt2 = 1, 4
	         do gen1 = 1, 3
	         do gen2 = 1, 3

          Aij(sfmeigenout1,sfmeigenin1,tt1,tt2,gen1,gen2) = (USf(sfmeigenout1,1,tt1,gen1)*USf(sfmeigenin1,1,tt1,gen1) - USf(sfmeigenout1,2,tt1,gen1)*USf(sfmeigenin1,2,tt1,gen1))*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)
		  enddo
		  enddo
		  enddo
	      enddo
	   enddo
     	enddo
! 	write(*,*)'------------------------------------------------------'
!        	write(*,*)'VSfSf(k,isf1,isf2,itt1,itt2,igen1,igen2)',VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)
! 	write(*,*)'ALt1(n,isf1,itt1,ftt1,igen1,fgen1)',ALt1(1,isf1,itt1,ftt1,igen1,fgen1)
! 	write(*,*)'ARt1(n,isf1,itt1,ftt1,igen1,fgen1)',ARt1(1,isf1,itt1,ftt1,igen1,fgen1)
! 	write(*,*)'ALt2(n,isf2,itt2,ftt2,igen2,fgen2)',ALt2(1,isf2,itt2,ftt2,igen2,fgen2)
! 	write(*,*)'ARt2(n,isf2,itt2,ftt2,igen2,fgen2)',ARt2(1,isf2,itt2,ftt2,igen2,fgen2)
! 	write(*,*)'BR(k,ftt1,ftt2,fgen1,fgen2)',BR(3,ftt1,ftt2,fgen1,fgen2)
! 	write(*,*)'BL(k,ftt1,ftt2,fgen1,fgen2)',BL(3,ftt1,ftt2,fgen1,fgen2)
! 	write(*,*)'------------------------------------------------------'
	end
