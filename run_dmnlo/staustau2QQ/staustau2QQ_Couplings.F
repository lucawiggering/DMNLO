      subroutine staustau2QQ_SetCouplings
      implicit none

#include "staustau2QQ_Model.h"
#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_GenCouplings.h"


      integer i, j, k, l, m, n, xtype, ytype, a, b, o, p, tt1, tt2, gen1, gen2,sfmeigenout1,sfmeigenin1

! functions from DMNLO_Kinematics.F

      double precision Krondelta, ef, I3f, hf


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
	
 
	 
!******************************************************** 
! Vertex Gluino  - Squark (outgoing) - Quark (incoming) *
!********************************************************
!  e.g. GluL/R(sfmeigenin,sftypein,ftypeout,sfgenin,fgenout), couplings multiplied by -i

           do a=1, 2                                                                               
             do b= 1, 4
               do m= 1, 4
                 do i= 1, 3
                    do n= 1, 3

                     if((b.gt.2).and.(m.gt.2)) then

                        GluL(a,b,m,i,n) = - GS*sqrt2*USf(a,1,b,i)*Krondelta(b,m)*Krondelta(i,n)
                        GluR(a,b,m,i,n) =   GS*sqrt2*USf(a,2,b,i)*Krondelta(b,m)*Krondelta(i,n)
                        GluLt2(a,b,m,i,n) = GluL(a,b,m,i,n)
                        GluRt2(a,b,m,i,n) = GluR(a,b,m,i,n)
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
          
!**************************	
! Fermion-Higgs couplings *
!**************************      
!e.g. DL/R(i,f3type,f4type,f3gen,f4gen), couplings multiplied by -i

        do i = 1, 4
           do j = 1, 4
              do m = 1, 3
                 do n = 1, 3
!                  DL(1,i,j,m,n) = 1d0
!                  DR(1,i,j,m,n) = 1d0
!            
!--- incoming up-type fermion ---           

	
           
              if ((i.eq.1).or.(i.eq.3)) then
! h              
                  DL(1,i,j,m,n) = -(EL*Mf(i,m)/(sqrt2*MW*SW*SB))*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(1,i,j,m,n) =  DL(1,i,j,m,n)
! H                  
                  DL(2,i,j,m,n) = -(EL*Mf(i,m)/(sqrt2*MW*SW*SB))*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(2,i,j,m,n) =  DL(2,i,j,m,n)

 
!                   
! --- incoming down-type fermion --- 

 	             else
! h 	
 	          DL(1,i,j,m,n) = (EL*Mf(i,m)/(sqrt2*MW*SW*SB))*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)   
 	          DR(1,i,j,m,n) = DL(1,i,j,m,n)
 	          
!  	          hff(1)= (-SA*EL*Mf(4,3))/(2*MW*SB*SW) 
! H 	          
 	          DL(2,i,j,m,n) = -(EL*Mf(i,m)/(sqrt2*MW*SW*SB))*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
 	          DR(2,i,j,m,n) =  DL(2,i,j,m,n)

           	          
      endif

      enddo
      enddo
      enddo
      enddo          
          
!**************************
! Higgs-Sfermion-Sfermion *
!**************************
!! couplings multiplied by -i    

!From Higgs-Hunter (p.394 f)

        do i= 1, 4
           do j=1, 3

!GLR(higgs,sfermionnrin,out,sfermiontype,sfgen)

	if ((i.eq.1).or.(i.eq.3)) then


	GLR(1,1,1,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*CA  + EL*MZ*(I3f(i)-ef(i)*SW2)*SAB/(SW*CW)               
	GLR(1,2,2,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*CA + EL*MZ*ef(i)*SW2*SAB/(SW*CW)                      
	GLR(1,2,1,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*CA + MUE*SA) 
	GLR(1,1,2,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*CA + MUE*SA) 

	GLR(2,1,1,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*SA - EL*MZ*(I3f(i)-ef(i)*SW2)*CAB/(SW*CW)             
	GLR(2,2,2,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*SA - EL*MZ*ef(i)*SW2*CAB/(SW*CW)                       
	GLR(2,2,1,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*SA - MUE*CA)                                           
	GLR(2,1,2,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*SA - MUE*CA)    
      
	GLR(3,1,1,i,j) = (0d0,0d0)
	GLR(3,2,2,i,j) = (0d0,0d0)
	GLR(3,2,1,i,j) = cI*hf(i,j)/sqrt2*(Af(i,j)*CB + MUE*SB) 
	GLR(3,1,2,i,j) = -cI*hf(i,j)/sqrt2*(Af(i,j)*CB + MUE*SB) 

	GLR(4,1,1,i,j) = (0d0,0d0)
	GLR(4,2,2,i,j) = (0d0,0d0)
	GLR(4,2,1,i,j) = cI*hf(i,j)/sqrt2*(Af(i,j)*SB - MUE*CB) 
	GLR(4,1,2,i,j) = -cI*hf(i,j)/sqrt2*(Af(i,j)*SB - MUE*CB) 

! G_{LR}^{tb} where i = up-type
	GLR(5,1,1,i,j) = hf(i+1,j)*Mf(i+1,j)*SB + hf(i,j)*Mf(i,j)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,i,j) = hf(i,j)*Mf(i+1,j)*CB + hf(i+1,j)*Mf(i,j)*SB 
	GLR(5,2,1,i,j) = hf(i,j)*(Af(i,j)*CB + MUE*SB)
	GLR(5,1,2,i,j) = hf(i+1,j)*(Af(i+1,j)*SB + MUE*CB)

	GLR(6,1,1,i,j) = -hf(i+1,j)*Mf(i+1,j)*CB + hf(i,j)*Mf(i,j)*SB + EL*MW/(sqrt2*SW)*C2B
	GLR(6,2,2,i,j) = hf(i,j)*Mf(i+1,j)*SB - hf(i+1,j)*Mf(i,j)*CB                                  
	GLR(6,2,1,i,j) = hf(i,j)*(Af(i,j)*SB - MUE*CB)
	GLR(6,1,2,i,j) = hf(i+1,j)*(-Af(i+1,j)*CB + MUE*SB)

	else
	

	GLR(1,1,1,i,j) = sqrt2*(EL*Mf(i,m)/(sqrt2*MW*SW*SB))*Mf(i,j)*SA + EL*MZ*(I3f(i)-ef(i)*SW2)*SAB/(SW*CW)
	GLR(1,2,2,i,j) = sqrt2*hf(i,j)*Mf(i,j)*SA + EL*MZ*ef(i)*SW2*SAB/(SW*CW) 
	GLR(1,2,1,i,j) = hf(i,j)/sqrt2*(Af(i,j)*SA + MUE*CA) 
	GLR(1,1,2,i,j) = hf(i,j)/sqrt2*(Af(i,j)*SA + MUE*CA) 

	GLR(2,1,1,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*CA - EL*MZ*(I3f(i)-ef(i)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,i,j) = -sqrt2*hf(i,j)*Mf(i,j)*CA - EL*MZ*ef(i)*SW2*CAB/(SW*CW) 
	GLR(2,2,1,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*CA - MUE*SA) 
	GLR(2,1,2,i,j) = -hf(i,j)/sqrt2*(Af(i,j)*CA - MUE*SA) 

	GLR(3,1,1,i,j) = (0d0,0d0)
	GLR(3,2,2,i,j) = (0d0,0d0)
	GLR(3,2,1,i,j) = cI*hf(i,j)/sqrt2*(Af(i,j)*SB + MUE*CB) 
	GLR(3,1,2,i,j) = -cI*hf(i,j)/sqrt2*(Af(i,j)*SB + MUE*CB) 

	GLR(4,1,1,i,j) = (0d0,0d0)
	GLR(4,2,2,i,j) = (0d0,0d0)
	GLR(4,2,1,i,j) = -cI*hf(i,j)/sqrt2*(Af(i,j)*CB - MUE*SB) 
	GLR(4,1,2,i,j) = cI*hf(i,j)/sqrt2*(Af(i,j)*CB - MUE*SB) 

! G_{LR}^{bt} where i = down-type
      GLR(5,1,1,i,j) = hf(i,j)*Mf(i,j)*SB + hf(i-1,j)*Mf(i-1,j)*CB - EL*MW/(sqrt2*SW)*S2B
      GLR(5,2,2,i,j) = hf(i-1,j)*Mf(i,j)*CB + hf(i,j)*Mf(i-1,j)*SB
      GLR(5,2,1,i,j) = hf(i,j)*(Af(i,j)*SB + MUE*CB)
            GLR(5,1,2,i,j) = hf(i-1,j)*(Af(i-1,j)*CB + MUE*SB) 

      GLR(6,1,1,i,j) = -hf(i,j)*Mf(i,j)*CB + hf(i-1,j)*Mf(i-1,j)*SB + EL*MW/(sqrt2*SW)*C2B
      GLR(6,2,2,i,j) = hf(i-1,j)*Mf(i,j)*SB - hf(i,j)*Mf(i-1,j)*CB                                 
      GLR(6,2,1,i,j) = hf(i,j)*(-Af(i,j)*CB + MUE*SB)
      GLR(6,1,2,i,j) = hf(i-1,j)*(Af(i-1,j)*SB - MUE*CB) 

      endif

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
	  
!!BOTTOM COUPLING

!            hff(1)= (-SA*EL*2.750102)/(2*MW*CB*SW) !!DIFFERENT MASS FROM MICROMEGAS
!            hff(2)= (-EL*2.408967*CA)/(2*MW*CB*SW)
           
!!TOP COUPLING           
           
      hff(1)=(-ca*EL*Mf(3,3))/(2*MW*sb*SW)
      hff(2)=(-EL*Mf(3,3)*sa)/(2*MW*sb*SW)

           hSf(1)= (EL*(2*(CB*MW**2*SAB*SW**2*USf(1,1,2,3)**2-CB*MW**2*SAB*SW**2+CW**2*Mf(2,3)**2*SA) - C2W*CB*MW**2*SAB*USf(1,1,2,3)**2-
     -      2*(Af(2,3)*SA*Mf(2,3)+CA*MUE*Mf(2,3))*CW**2*USf(1,1,2,3)*USf(2,1,2,3)))/(2*CB*CW**2*MW*SW)
           hSf(2)= (EL*(-(2*(CA*CW**2*Mf(2,3)**2+CAB*CB*MW**2*SW**2*USf(1,1,2,3)**2-CAB*CB*MW**2*SW**2)-C2W*CAB*CB*MW**2*USf(1,1,2,3)**2-
     -         2*(Af(2,3)*CA*Mf(2,3)-MUE*SA*Mf(2,3))*CW**2*USf(1,1,2,3)*USf(2,1,2,3))))/(2*CB*CW**2*MW*SW)
            
             

	end
