	subroutine QQ2xx_SetCouplings
	implicit none

#include "QQ2xx_Model.h"
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"

	double complex ZLR(2,2,4,3),C(2,4),Xil,Xkj,Xij,Xkl 
	integer i,j,k,l,a,b,m,n,o,p,asftype1,asftype2,asfgen1,asfmeigen1,asfmeigen2,axhtype,ayhtype  
	double precision test
	
! function
	double precision Krondelta

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
                                
                BL(2,i,j,m,n) = -ef_(i)*EL*Krondelta(i,j)*Krondelta(m,n)
                BR(2,i,j,m,n) = -ef_(i)*EL*Krondelta(i,j)*Krondelta(m,n)                

                BV(i) = I3f_(i) - 2d0*SW**2*ef_(i)
                BA(i) = I3f_(i)

                BL(3,i,j,m,n) = - (EL/(2d0*CW*SW))*(BV(i) + BA(i))*Krondelta(i,j)*Krondelta(m,n)
                BR(3,i,j,m,n) = - (EL/(2d0*CW*SW))*(BV(i) - BA(i))*Krondelta(i,j)*Krondelta(m,n)

                if((i.eq.1).or.(i.eq.3)) then               !incoming up type (->W+ -radiation)

                BL(4,i,j,m,n) = - EL/(sqrt2*SW)*Krondelta(i,j-1)*Krondelta(m,n)

                else                                        !incoming down type (->W- -radiation)

                BL(4,i,j,m,n) = - EL/(sqrt2*SW)*Krondelta(i,j+1)*Krondelta(m,n)

                endif

                BR(4,i,j,m,n) = (0d0,0d0)

                do k=1,4  

                    if((k.gt.3).and.((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then
                      
                       if((i.eq.1).or.(i.eq.3)) then

                       BLt(k,i,j,m,n) = BL(k,i,j,m,n)
                       BLu(k,i,j,m,n) = 0d0

                       BRt(k,i,j,m,n) = BR(k,i,j,m,n)
                       BRu(k,i,j,m,n) = 0d0

                       else

                       BLu(k,i,j,m,n) = BL(k,i,j,m,n)
                       BLt(k,i,j,m,n) = 0d0

                       BRu(k,i,j,m,n) = BR(k,i,j,m,n)
                       BRt(k,i,j,m,n) = 0d0

                       endif

                    else

                       BLt(k,i,j,m,n) = BL(k,i,j,m,n)
                       BLu(k,i,j,m,n) = BL(k,i,j,m,n)

                       BRt(k,i,j,m,n) = BR(k,i,j,m,n)
                       BRu(k,i,j,m,n) = BR(k,i,j,m,n)  

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
           
! --- incoming up-type fermion ---           
           
              if ((i.eq.1).or.(i.eq.3)) then
! h              
                  DL(1,i,j,m,n) = -hf_(i,m)*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(1,i,j,m,n) =  DL(1,i,j,m,n)
! H                  
                  DL(2,i,j,m,n) = -hf_(i,m)*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(2,i,j,m,n) =  DL(2,i,j,m,n)
! H3                  
                  DL(3,i,j,m,n) = -cI*hf_(i,m)*CB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
                  DR(3,i,j,m,n) = -DL(3,i,j,m,n)
! G0                  
                  DL(4,i,j,m,n) = -cI*hf_(i,m)*SB/sqrt2*Krondelta(i,j)*Krondelta(m,n)
                  DR(4,i,j,m,n) = -DL(4,i,j,m,n)
! H+(out)                
		  DL(5,i,j,m,n) =  hf_(i+1,m)*SB*Krondelta(i,j-1)*Krondelta(m,n)                  
		  DR(5,i,j,m,n) =  hf_(i,m)*CB*Krondelta(i,j-1)*Krondelta(m,n)
! G+(out)                 
	          DL(6,i,j,m,n) = -hf_(i+1,m)*CB*Krondelta(i,j-1)*Krondelta(m,n)
		  DR(6,i,j,m,n) = hf_(i,m)*SB*Krondelta(i,j-1)*Krondelta(m,n)
 
                  do k =1,6
                    
                     if((k.gt.4).and.((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then
                   
                       DLt(k,i,j,m,n) = DL(k,i,j,m,n)
                       DLu(k,i,j,m,n) = 0d0

                       DRt(k,i,j,m,n) = DR(k,i,j,m,n)
                       DRu(k,i,j,m,n) = 0d0

                     else

                       DLt(k,i,j,m,n) = DL(k,i,j,m,n)
                       DLu(k,i,j,m,n) = DL(k,i,j,m,n)

                       DRt(k,i,j,m,n) = DR(k,i,j,m,n)
                       DRu(k,i,j,m,n) = DR(k,i,j,m,n)        

                     endif
     
   
                  enddo

! --- incoming down-type fermion --- 

 	             else
! h 	
 	          DL(1,i,j,m,n) = hf_(i,m)*SA/sqrt2*Krondelta(i,j)*Krondelta(m,n)   
 	          DR(1,i,j,m,n) = DL(1,i,j,m,n)
! H 	          
 	          DL(2,i,j,m,n) = -hf_(i,m)*CA/sqrt2*Krondelta(i,j)*Krondelta(m,n)
 	          DR(2,i,j,m,n) =  DL(2,i,j,m,n)
! H3 	          
 	          DL(3,i,j,m,n) = -cI*hf_(i,m)*SB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
 	          DR(3,i,j,m,n) = -DL(3,i,j,m,n)
! G0 	          
 	          DL(4,i,j,m,n) = cI*hf_(i,m)*CB/sqrt2*Krondelta(i,j)*Krondelta(m,n) 
 	          DR(4,i,j,m,n) = -DL(4,i,j,m,n)
! H-(out) 	          
 	          DL(5,i,j,m,n) = hf_(i-1,m)*CB*Krondelta(i,j+1)*Krondelta(m,n)
		  DR(5,i,j,m,n) = hf_(i,m)*SB*Krondelta(i,j+1)*Krondelta(m,n) 
! G-(out) 	          
		  DL(6,i,j,m,n) = hf_(i-1,m)*SB*Krondelta(i,j+1)*Krondelta(m,n)                  
		  DR(6,i,j,m,n) = -hf_(i,m)*CB*Krondelta(i,j+1)*Krondelta(m,n)
 
                  do k =1,6


                     if((k.gt.4).and.((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then
                   
                       DLu(k,i,j,m,n) = DL(k,i,j,m,n)
                       DLt(k,i,j,m,n) = 0d0

                       DRu(k,i,j,m,n) = DR(k,i,j,m,n)
                       DRt(k,i,j,m,n) = 0d0

                     else

                       DLt(k,i,j,m,n) = DL(k,i,j,m,n)
                       DLu(k,i,j,m,n) = DL(k,i,j,m,n)

                       DRt(k,i,j,m,n) = DR(k,i,j,m,n)
                       DRu(k,i,j,m,n) = DR(k,i,j,m,n)        

                     endif
   
                  enddo
 
 	           	          
 	          endif

                enddo
             enddo
 	   enddo
 	enddo

!******************
!Higgs-Higgs-Higgs*               
!******************
!e.g I(i,xhtype,yhtype), couplings multiplied by -i

!convention: 1=h, 2=H, 3=P, 4=G0, 5=H+, 6=G+

      do i=1, 6
        do j=1, 6
          do l=1, 6

           ILR(i,j,l) = (0d0,0d0)

          enddo
        enddo
      enddo



!h-h-h

      ILR(1,1,1)=-3*EL*MZ/(2*CW*SW)*C2A*SAB

!h-h-H

      ILR(1,1,2)=-EL*MZ/(2*CW*SW)*(2*S2A*SAB-CAB*C2A)

      ILR(1,2,1)=ILR(1,1,2)                            
      ILR(2,1,1)=ILR(1,1,2)

!h-H-H

      ILR(1,2,2)= EL*MZ/(2*CW*SW)*(2*S2A*CAB+SAB*C2A)

      ILR(2,2,1)=ILR(1,2,2)                            
      ILR(2,1,2)=ILR(1,2,2)

!h-P-P

      ILR(1,3,3)=-EL*MZ/(2*CW*SW)*C2B*SAB

      ILR(3,3,1)=ILR(1,3,3)                           
      ILR(3,1,3)=ILR(1,3,3)

!h-P-G0

      ILR(1,3,4)=-EL*MZ/(2*SW*CW)*S2B*SAB

      ILR(1,4,3)=ILR(1,3,4)                           
      ILR(3,1,4)=ILR(1,3,4)
      ILR(3,4,1)=ILR(1,3,4)
      ILR(4,1,3)=ILR(1,3,4)
      ILR(4,3,1)=ILR(1,3,4)
     
!h-G0-G0

      ILR(1,4,4)=EL*MZ/(2*CW*SW)*C2B*SAB              

      ILR(4,4,1)=ILR(1,4,4)                            
      ILR(4,1,4)=ILR(1,4,4)

!h-H+-H-

      ILR(1,5,5)=-EL/SW*(MW*SBA+MZ/(2*CW)*C2B*SAB)

      ILR(5,5,1)=ILR(1,5,5)                            
      ILR(5,1,5)=ILR(1,5,5)

!h-G+-H-

      ILR(1,5,6)=EL/(2*SW)*(MW*CBA-MZ/CW*S2B*SAB)

      ILR(1,6,5)=ILR(1,5,6)                            
      ILR(5,1,6)=ILR(1,5,6)
      ILR(5,6,1)=ILR(1,5,6)
      ILR(6,1,5)=ILR(1,5,6)
      ILR(6,5,1)=ILR(1,5,6)

!h-G+-G-

      ILR(1,6,6)=EL/(2*CW*SW)*MZ*C2B*SAB             

      ILR(6,6,1)=ILR(1,6,6)                            
      ILR(6,1,6)=ILR(1,6,6)

!H-H-H

      ILR(2,2,2)=-3*EL*MZ/(2*CW*SW)*C2A*CAB

!H-P-P

      ILR(2,3,3)=EL*MZ/(2*CW*SW)*C2B*CAB

      ILR(3,3,2)=ILR(2,3,3)                            
      ILR(3,2,3)=ILR(2,3,3)

!H-P-G0

      ILR(2,3,4)=EL*MZ/(2*SW*CW)*S2B*CAB

      ILR(2,4,3)=ILR(2,3,4)                            
      ILR(3,2,4)=ILR(2,3,4)
      ILR(3,4,2)=ILR(2,3,4)
      ILR(4,2,3)=ILR(2,3,4)
      ILR(4,3,2)=ILR(2,3,4)

!H-G0-G0

      ILR(2,4,4)=-EL*MZ/(2*CW*SW)*C2B*CAB              

      ILR(4,4,2)=ILR(2,4,4)                            
      ILR(4,2,4)=ILR(2,4,4)

!H-H+-H-

      ILR(2,5,5)=-EL/SW*(MW*CBA-MZ/(2*CW)*C2B*CAB)

      ILR(5,5,2)=ILR(2,5,5)                            
      ILR(5,2,5)=ILR(2,5,5)

!H-G+-H-

      ILR(2,5,6)=-EL/(2*SW)*(MW*SBA-MZ/CW*S2B*CAB)

      ILR(2,6,5)=ILR(2,5,6)                            
      ILR(5,2,6)=ILR(2,5,6)
      ILR(5,6,2)=ILR(2,5,6)
      ILR(6,2,5)=ILR(2,5,6)
      ILR(6,5,2)=ILR(2,5,6)

!H-G+-G-

      ILR(2,6,6)=-EL/(2*CW*SW)*MZ*C2B*CAB     

      ILR(6,2,6)=ILR(2,6,6)                            
      ILR(6,6,2)=ILR(2,6,6)

!P-H+(in)-G+(out)                                            
!always used as ILR(channel,xtype,ytype)

      ILR(3,5,6)=cI*EL/(2*SW)*MW                       !P(channel),H+(out),G-(out)
      ILR(3,6,5)=-ILR(3,5,6)                            
      ILR(5,3,6)=ILR(3,5,6)
      ILR(5,6,3)=-ILR(3,5,6)
      ILR(6,3,5)=-ILR(3,5,6)
      ILR(6,5,3)=ILR(3,5,6)

!********************************
! Vector-Vector-Higgs couplings *
!********************************
!e.g. LLR(2,2,yhtype) = Z-Z-Higgs, couplings multiplied by -i

!1=Gluon (for convention) 2=Gamma, 3=Z-Boson, 4=W-Boson

      do i=1, 4
         do j=1, 4
            do k= 1, 6

      LLR(i,j,k) = (0d0,0d0)
     
            enddo
         enddo
      enddo

      
!W+-W--H

      LLR(4,4,2)=EL*MZ/TW*CBA

!W+-W--h

      LLR(4,4,1)=EL*MZ/TW*SBA

!Z-Z-H

      LLR(3,3,2)=2*EL*MZ/S2W*CBA

!Z-Z-h

      LLR(3,3,1)=2*EL*MZ/S2W*SBA

!G-W+-G+

      LLR(2,4,6)=EL*MW
      LLR(4,2,6)=LLR(2,4,6)

!Z-W+-G+

      LLR(3,4,6)=-EL*MZ*SW      
      LLR(4,3,6)=LLR(3,4,6)

!******************************
! Ghost-Ghost-Higgs couplings *
!******************************
!e.g.GGH(2,2,yhtype) = G_Z-G_Z-Higgs, couplings multiplied by -i

!1=Gluon (for convention) 2=Gamma, 3=Z-Boson, 4=W-Boson  , from Higgs-Hunter p. 385 and CalcHEP-MSSM-model-file + Rosiek p. 36 + Feynarts-MSSM-model p.39

      do i=1, 4
         do j=1, 4
            do k= 1, 6   

      GGH(i,j,k) = (0d0,0d0)
     
            enddo
         enddo
      enddo

!Gh_Z-Gh_Z-h

      GGH(3,3,1)=-EL*MZ/S2W*SBA

!Gh_Z-Gh_Z-H

      GGH(3,3,2)=-EL*MZ/S2W*CBA

!Gh_W+-Gh_W--h

      GGH(4,4,1)=-EL*MZ/(2d0*TW)*SBA   

!Gh_W+-Gh_W--H

      GGH(4,4,2)=-EL*MZ/(2d0*TW)*CBA   

!Gh_W+-Gh_W--G0

      GGH(4,4,4)=-cI*EL*MZ/(2d0*TW)    


      if(total_chargeout.gt.0) then

!Gh_W+(out)-Gh_G(out)-G+(in)

      GGH(4,2,6)=-EL*MW            

!Gh_W+(in)-Gh_G(in)-G+(out)           !coupling to use instead of the hermitian conjugated of GGH(4,2,6)

      GGH(2,4,6)= 0d0                 !coupling doesn't exist         

!Gh_W+(out)-Gh_Z(out)-G+(in)  

      GGH(4,3,6)=-EL*MW*C2W/S2W

!Gh_W+(in)-Gh_Z(in)-G+(out)           !coupling to use instead of the hermitian conjugated of GGH(4,3,6)

      GGH(3,4,6)= EL*MW/(S2W) 

      else    

!Gh_W-(in)-Gh_G(in)-G-(out) 

      GGH(4,2,6)= 0d0                 !coupling doesn't exist

!Gh_W-(out)-Gh_G(out)-G-(in)          !coupling to use instead of the hermitian conjugated of GGH(4,2,6)

      GGH(2,4,6)=-EL*MW               !typo in Rosiek p. 37 (wrong sign,see for example Lagrangian in Rosiek p.15  + Feynarts-MSSM-model p.39 for correct version)     

!Gh_W-(in)-Gh_Z(in)-G-(out)

      GGH(4,3,6)= EL*MW/(S2W)

!Gh_W-(out)-Gh_Z(out)-G-(in)          !coupling to use instead of the hermitian conjugated  of GGH(4,3,6) 

      GGH(3,4,6)=-EL*MW*C2W/S2W    

      endif

!*******************************
! Ghost-Ghost-Vector couplings *
!*******************************
!e.g.GGV(2,2,yhtype) = G_Z-G_Z-Higgs, couplings multiplied by -i

!1=Gluon (for convention) 2=Gamma, 3=Z-Boson, 4=W-Boson, taken from CalcHEP-MSSM-model-file + Rosiek (p.36 + p.46)

      do i=1, 4
         do j=1, 4
            do k= 1, 4   

      GGV(i,j,k) = (0d0,0d0)
     
            enddo
         enddo
      enddo


!Gh_G-Gh_G--G                              !Peskin p. 515+ eq.16.42 + fig.16.2      

      GGV(1,1,1) = -cI*GS                  !necessary sign-change (see e.g. Rosiek p.46) included via using the hermitian conjugated couplings in conjugated diagramms (all other GGV-couplings are real)

!Gh_W+(out)-Gh_W+(in)-g

      GGV(4,4,2) = EL

!Gh_W+(out)-Gh_W+(in)-Z

      GGV(4,4,3) = EL/TW

!Gh_W+(out)-Gh_G(out)-W+(in)

      GGV(4,2,4) = -EL
      GGV(2,4,4) = -GGV(4,2,4)             !needs to be used instead of the hermitian conjugated coupling in conjugated diagramms

!Gh_W+(out)-Gh_Z(out)-W+(in)

      GGV(4,3,4) = -EL/TW
      GGV(3,4,4) = -GGV(4,3,4)             !needs to be used instead of the hermitian conjugated coupling in conjugated diagramms


!*********************************
! Vector-Vector-Vector couplings *
!*********************************
!e.g. VVV(i,xtype,ytype), couplings multiplied by -i

!1=Gluon 2=Gamma, 3=Z-Boson, 4=W+(W-)-Boson, couplings taken from Sparticles and FeynArts MSSM.dvi and MSSMQCD.dvi + for conventions used in kinematics see Mathematics-file

      do i=1, 4
         do j=1, 4
            do k= 1, 4

      VVV(i,j,k) = (0d0,0d0)
     
            enddo
         enddo
      enddo

!G-G-G                                                   

      VVV(1,1,1) = cI*GS       

!g(in)-W+(in)-W-(in)

      VVV(2,4,4) = -EL
      VVV(4,2,4) =-VVV(2,4,4)
      VVV(4,4,2) =-VVV(2,4,4)


!Z(in)-W+(in)-W-(in)

      VVV(3,4,4) = -EL/TW
      VVV(4,3,4) =-VVV(3,4,4)
      VVV(4,4,3) =-VVV(3,4,4)


!*******************************
! Vector-Higgs-Higgs couplings *
!*******************************
!e.g. JLR(2,xhtype(in),yhtype(out)) = Z-h-P = -JLR(2,xhtype(out),yhtype(in)), couplings multiplied by -i

!1=Gluon 2=Gamma, 3=Z-Boson, 4=W-Boson

!conventions used in kinematics: HH~JLR*(k1-k2), VH~JLR*(p1+p2+k2), HV~-JLR*(p1+p2+k1)

      do i=1, 4
         do j=1, 6
            do k= 1, 6

      JLR(i,j,k) = (0d0,0d0)
     
            enddo
         enddo
      enddo 

!G-H-(out)-H-(in)    
 
      if(total_chargeout.eq.0) then      

      JLR(2,5,5)=-EL

      else

      JLR(2,5,5)=EL

      endif

!G-G-(out)-G-(in)

      if(total_chargeout.eq.0) then

      JLR(2,6,6)=-EL    

      else  

      JLR(2,6,6)=EL 

      endif              

!Z-h(in)-P(out)

      JLR(3,1,3)=-cI*EL/(2*CW*SW)*CBA      

      JLR(3,3,1)=-JLR(3,1,3)               !sign change from momenta/kinematics for h(out),P(in)

!Z-h(in)-G0(out)

      JLR(3,1,4)=-cI*EL/(2*CW*SW)*SBA    

      JLR(3,4,1)=-JLR(3,1,4)               !sign change from momenta/kinematics for h(out),G0(in)

!Z-H(in)-P(out)

      JLR(3,2,3)=cI*EL/(2*CW*SW)*SBA

      JLR(3,3,2)=-JLR(3,2,3)               !sign change from momenta/kinematics for H(out),P(in)

!Z-H(in)-G0(out)

      JLR(3,2,4)=-cI*EL/(2*CW*SW)*CBA     

      JLR(3,4,2)=-JLR(3,2,4)               !sign change from momenta/kinematics for H(out),G0(in)

!Z-H-(out)-H-(in)=-(Z-H+(out)-H+(in))

      if(total_chargeout.eq.0) then

      JLR(3,5,5)=-EL*C2W/(2*SW*CW)

      else

      JLR(3,5,5)=EL*C2W/(2*SW*CW)

      endif

!Z-G-(out)-G-(in)

      if(total_chargeout.eq.0) then

      JLR(3,6,6)=-EL*C2W/(2*SW*CW)  

      else

      JLR(3,6,6)=EL*C2W/(2*SW*CW)

      endif
     

!W-(in)-h-H-(out)                          !kinematics different for in/out-choice (Sparticles S. 523, higgs-Hunter S. 358)

      if(total_chargeout.eq.0) then

      JLR(4,1,5)= EL/(2*SW)*CBA          

      JLR(4,5,1)= JLR(4,1,5)   

      else

      JLR(4,1,5)=-EL/(2*SW)*CBA          

      JLR(4,5,1)= JLR(4,1,5) 

      endif          

!W-(in)-h-G-(out)                          !kinematics different for in/out-choice (Sparticles S 523, higgs-Hunter S. 358)  

      if(total_chargeout.eq.0) then

      JLR(4,1,6)= EL/(2*SW)*SBA

      JLR(4,6,1)= JLR(4,1,6)

      else

      JLR(4,1,6)=-EL/(2*SW)*SBA

      JLR(4,6,1)= JLR(4,1,6)

      endif
 
!W-(in)-H-H-(out)                          !kinematics different for in/out-choice (Sparticles S 523, higgs-Hunter S. 358)

      if(total_chargeout.eq.0) then

      JLR(4,2,5)=-EL/(2*SW)*SBA        
 
      JLR(4,5,2)= JLR(4,2,5)

      else

      JLR(4,2,5)= EL/(2*SW)*SBA        
 
      JLR(4,5,2)= JLR(4,2,5) 

      endif 

!W-(in)-H-G-(out)                          !kinematics different for in/out-choice (Sparticles S 523, higgs-Hunter S. 358)

      if(total_chargeout.eq.0) then

      JLR(4,2,6)= EL/(2*SW)*CBA       
 
      JLR(4,6,2)= JLR(4,2,6)  

      else  

      JLR(4,2,6)=-EL/(2*SW)*CBA       
 
      JLR(4,6,2)= JLR(4,2,6) 

      endif      
    
!W+-(in)-P(in)-H+-(out)                    !kinematics different for in/out-choice (Sparticles S 523, higgs-Hunter S. 358)

       JLR(4,3,5)=-cI*EL/(2*SW)          

       JLR(4,5,3)=-JLR(4,3,5)              !sign change from momenta/kinematics for G+(in),H(out) 

!W+-(in)-G0-G+-(out)                       !kinematics different for different in/out-choice (Sparticles S 523, higgs-Hunter S. 358)

       JLR(4,4,6)=-cI*EL/(2*SW)         

       JLR(4,6,4)=-JLR(4,4,6)              !sign change from momenta/kinematics for G+(in),H(out) 


!******************************************
! sFermion-sFermion-Higgs-Higgs-couplings *
!******************************************
!e.g. HHSfSf(xhtype,yhtype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1), no flavor mixing => sfgen2-independent, couplings multiplied by -i

!neglect flavor-mixing
!used Higgs-Hunter

!1=h, 2=H, 3=P, 4= G0, 5=Hp, 6= Gp


        do axhtype=1,6                 !no flavor-mixing, hence sfgen2 not needed
          do ayhtype=1,6
            do asfmeigen1=1,2
              do asfmeigen2= 1,2
                do asftype1=1,4
                  do asfgen1=1,3

                   KLR(axhtype,ayhtype,asfmeigen1,asfmeigen2,asftype1,asfgen1)= (0d0,0d0)

                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do i= 1, 4                    
           do j=1, 3

	if ((i.eq.1).or.(i.eq.3)) then
       
!h,h

	KLR(1,1,1,1,i,j) = -hf_(i,j)**2*CA2 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2A                      !see Higgs-Hunter p. 404
	KLR(1,1,2,2,i,j) = -hf_(i,j)**2*CA2 + EL**2/(2d0*CW2)*ef_(i)*C2A                    

!H,H
 
        KLR(2,2,1,1,i,j) = -hf_(i,j)**2*SA2 - 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2A    
	KLR(2,2,2,2,i,j) = -hf_(i,j)**2*SA2 - EL**2/(2d0*CW2)*ef_(i)*C2A    

!P,P

        KLR(3,3,1,1,i,j) = -hf_(i,j)**2*CB2 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2B    
	KLR(3,3,2,2,i,j) = -hf_(i,j)**2*CB2 + EL**2/(2d0*CW2)*ef_(i)*C2B           

!G0,G0

        KLR(4,4,1,1,i,j) = -hf_(i,j)**2*SB2 - 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2B    
	KLR(4,4,2,2,i,j) = -hf_(i,j)**2*SB2 - EL**2/(2d0*CW2)*ef_(i)*C2B     
             
!h,H

	KLR(1,2,1,1,i,j) = -hf_(i,j)**2*S2A/2d0 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*S2A                 
	KLR(1,2,2,2,i,j) = -hf_(i,j)**2*S2A/2d0 + EL**2/(2d0*CW2)*ef_(i)*S2A                               
	kLR(2,1,1,1,i,j) = KLR(1,2,1,1,i,j)                                                              
	KLR(2,1,2,2,i,j) = KLR(1,2,2,2,i,j)                                                                          

!G0,P
      
	KLR(3,4,1,1,i,j) = -hf_(i,j)**2*S2B/2d0 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*S2B 
	KLR(3,4,2,2,i,j) = -hf_(i,j)**2*S2B/2d0 + EL**2/(2d0*CW2)*ef_(i)*S2B
	KLR(4,3,1,1,i,j) = KLR(3,4,1,1,i,j)
	KLR(4,3,2,2,i,j) = KLR(3,4,2,2,i,j)

!Hp-Hm

	KLR(5,5,1,1,i,j) = EL**2/(2d0*SW2)*C2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2) - 2*EL**2/(MZ2*S2W**2)*Mf2(i+1,j)*TB2
	KLR(5,5,2,2,i,j) = EL**2/(2d0*CW2)*C2B*ef_(i)-2*EL**2*Mf2(i,j)/(MZ2*S2W**2*TB2)

!Gp-Gm

	KLR(6,6,1,1,i,j) = -EL**2/(2d0*SW2)*C2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2) - 2*EL**2/(MZ2*S2W**2)*Mf2(i+1,j)
	KLR(6,6,2,2,i,j) = -EL**2/(2d0*CW2)*C2B*ef_(i)-2*EL**2*Mf2(i,j)/(MZ2*S2W**2)

!Gp,Hm (Gm,Hp)

	KLR(6,5,1,1,i,j) = EL**2/(2d0*SW2)*S2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2) + 2*EL**2/(MZ2*S2W**2)*Mf2(i+1,j)*TB
	KLR(6,5,2,2,i,j) = EL**2/(2d0*SW2)*S2B*ef_(i)*TW2 - EL**2/(2d0*SW2*MW2*TB)*Mf2(i,j)
	KLR(5,6,1,1,i,j) = KLR(6,5,1,1,i,j)
	KLR(5,6,2,2,i,j) = KLR(6,5,2,2,i,j)

!neutral Higgs,Hp (Hm)

!error in the HiggsHunter p.408, necessary replacements: SA->CA,CA->-SA 
	KLR(1,5,1,1,i,j) =-EL**2/(2d0*sqrt2*SW2)*(CAB-Mf2(i,j)/MW2*CA*CB/SB2-Mf2(i+1,j)/MW2*(-SA)*SB/CB2)   
	KLR(1,5,2,2,i,j) = hf_(i,j)*hf_(i+1,j)/sqrt2*SBA
	KLR(5,1,1,1,i,j) = KLR(1,5,1,1,i,j)
	KLR(5,1,2,2,i,j) = KLR(1,5,2,2,i,j)

!error in the HiggsHunter p.408, necessary replacements: SA->-CA,CA->SA
	KLR(2,5,1,1,i,j) =-EL**2/(2d0*sqrt2*SW2)*(SAB-Mf2(i,j)/MW2*SA*CB/SB2+Mf2(i+1,j)/MW2*(-CA)*SB/CB2)
	KLR(2,5,2,2,i,j) = hf_(i,j)*hf_(i+1,j)/sqrt2*CBA
	KLR(5,2,1,1,i,j) = KLR(2,5,1,1,i,j)
	KLR(5,2,2,2,i,j) = KLR(2,5,2,2,i,j)

	KLR(3,5,1,1,i,j) =-cI*EL**2/(2d0*sqrt2*SW2)*(C2B-Mf2(i,j)/(MW2*TB2)+Mf2(i+1,j)/MW2*TB2) 
	KLR(3,5,2,2,i,j) = (0d0,0d0)
	KLR(5,3,1,1,i,j) = KLR(3,5,1,1,i,j)
	KLR(5,3,2,2,i,j) = KLR(3,5,2,2,i,j)

	KLR(4,5,1,1,i,j) =-cI*EL**2/(2d0*sqrt2*SW2)*(S2B-Mf2(i,j)/(MW2*TB)-Mf2(i+1,j)/MW2*TB) 
	KLR(4,5,2,2,i,j) =-cI*hf_(i,j)*hf_(i+1,j)/sqrt2
	KLR(5,4,1,1,i,j) = KLR(4,5,1,1,i,j)
	KLR(5,4,2,2,i,j) = KLR(4,5,2,2,i,j)

!neutral Higgs,Gp (Gm)

!error in the HiggsHunter p.408, necessary replacements: SA->-CA,CA->SA
	KLR(1,6,1,1,i,j) = EL**2/(2d0*sqrt2*SW2)*(-SAB-Mf2(i,j)/MW2*(-CA)/SB+Mf2(i+1,j)/MW2*SA/CB) 
	KLR(1,6,2,2,i,j) = -hf_(i,j)*hf_(i+1,j)/sqrt2*CBA
	KLR(6,1,1,1,i,j) = KLR(1,6,1,1,i,j)
	KLR(6,1,2,2,i,j) = KLR(1,6,2,2,i,j)

!error in the HiggsHunter p.408, necessary replacements: SA->CA,CA->-SA
	KLR(2,6,1,1,i,j) = EL**2/(2d0*sqrt2*SW2)*(CAB-Mf2(i,j)/MW2*(-SA)/SB-Mf2(i+1,j)/MW2*CA/CB)
	KLR(2,6,2,2,i,j) = hf_(i,j)*hf_(i+1,j)/sqrt2*SBA
	KLR(6,2,1,1,i,j) = KLR(2,6,1,1,i,j)
	KLR(6,2,2,2,i,j) = KLR(2,6,2,2,i,j)

	KLR(3,6,1,1,i,j) =-cI*EL**2/(2d0*sqrt2*SW2)*(S2B-Mf2(i,j)/(MW2*TB)-Mf2(i+1,j)/MW2*TB) 
	KLR(3,6,2,2,i,j) = cI*hf_(i,j)*hf_(i+1,j)/sqrt2
	KLR(6,3,1,1,i,j) = KLR(3,6,1,1,i,j)
	KLR(6,3,2,2,i,j) = KLR(3,6,2,2,i,j)

	KLR(4,6,1,1,i,j) = cI*EL**2/(2d0*sqrt2*SW2)*(C2B+Mf2(i,j)/MW2-Mf2(i+1,j)/MW2) 
	KLR(4,6,2,2,i,j) = (0d0,0d0)
	KLR(6,4,1,1,i,j) = KLR(4,6,1,1,i,j)
	KLR(6,4,2,2,i,j) = KLR(4,6,2,2,i,j)
	else

!h,h

	KLR(1,1,1,1,i,j) = -hf_(i,j)**2*SA2 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2A                      !see Higgs-Hunter p. 404
	KLR(1,1,2,2,i,j) = -hf_(i,j)**2*SA2 + EL**2/(2d0*CW2)*ef_(i)*C2A                    

!H,H
 
        KLR(2,2,1,1,i,j) = -hf_(i,j)**2*CA2 - 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2A    
	KLR(2,2,2,2,i,j) = -hf_(i,j)**2*CA2 - EL**2/(2d0*CW2)*ef_(i)*C2A    

!P,P

        KLR(3,3,1,1,i,j) = -hf_(i,j)**2*SB2 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2B    
	KLR(3,3,2,2,i,j) = -hf_(i,j)**2*SB2 + EL**2/(2d0*CW2)*ef_(i)*C2B           

!G0,G0

        KLR(4,4,1,1,i,j) = -hf_(i,j)**2*CB2 - 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*C2B    
	KLR(4,4,2,2,i,j) = -hf_(i,j)**2*CB2 - EL**2/(2d0*CW2)*ef_(i)*C2B        

!h,H

	KLR(1,2,1,1,i,j) = hf_(i,j)**2*S2A/2d0 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*S2A                   
	KLR(1,2,2,2,i,j) = hf_(i,j)**2*S2A/2d0 + EL**2/(2d0*CW2)*ef_(i)*S2A                                
	kLR(2,1,1,1,i,j) = KLR(1,2,1,1,i,j)                                                              
	KLR(2,1,2,2,i,j) = KLR(1,2,2,2,i,j)                                                                          

!G0,P
      
	KLR(3,4,1,1,i,j) = hf_(i,j)**2*S2B/2d0 + 2*EL**2/S2W**2*(I3f_(i)-ef_(i)*SW2)*S2B 
	KLR(3,4,2,2,i,j) = hf_(i,j)**2*S2B/2d0 + EL**2/(2d0*CW2)*ef_(i)*S2B
	KLR(4,3,1,1,i,j) = KLR(3,4,1,1,i,j)
	KLR(4,3,2,2,i,j) = KLR(3,4,2,2,i,j)

!Hp-Hm

	KLR(5,5,1,1,i,j) = EL**2/(2d0*SW2)*C2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2)-2*EL**2/(MZ2*S2W**2)*Mf2(i-1,j)/TB2
	KLR(5,5,2,2,i,j) = EL**2/(2d0*CW2)*C2B*ef_(i)-2*EL**2*Mf2(i,j)/(MZ2*S2W**2)*TB2

!Gp-Gm

	KLR(6,6,1,1,i,j) = -EL**2/(2d0*SW2)*C2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2)-2*EL**2/(MZ2*S2W**2)*Mf2(i-1,j)
	KLR(6,6,2,2,i,j) = -EL**2/(2d0*CW2)*C2B*ef_(i)-2*EL**2*Mf2(i,j)/(MZ2*S2W**2)

!Gp,Hm (Gm,Hp)

	KLR(6,5,1,1,i,j) = EL**2/(2d0*SW2)*S2B*(-2*I3f_(i)+(I3f_(i)-ef_(i)*SW2)/CW2) - 2*EL**2/(MZ2*S2W**2)*Mf2(i-1,j)/TB
	KLR(6,5,2,2,i,j) = EL**2/(2d0*SW2)*S2B*ef_(i)*TW2 + EL**2*TB/(2d0*SW2*MW2)*Mf2(i,j)
	KLR(5,6,1,1,i,j) = KLR(6,5,1,1,i,j)
	KLR(5,6,2,2,i,j) = KLR(6,5,2,2,i,j)

!neutral Higgs,Hp (Hm)

!error in the HiggsHunter p.408, necessary replacements: SA->CA,CA->-SA (wrong->correct)
	KLR(1,5,1,1,i,j) =-EL**2/(2d0*sqrt2*SW2)*(CAB-Mf2(i-1,j)/MW2*CA*CB/SB2-Mf2(i,j)/MW2*(-SA)*SB/CB2) 
	KLR(1,5,2,2,i,j) = hf_(i-1,j)*hf_(i,j)/sqrt2*SBA
	KLR(5,1,1,1,i,j) = KLR(1,5,1,1,i,j)
	KLR(5,1,2,2,i,j) = KLR(1,5,2,2,i,j)

!error in the HiggsHunter p.408, necessary replacements: SA->-CA,CA->SA (wrong->correct)
	KLR(2,5,1,1,i,j) =-EL**2/(2d0*sqrt2*SW2)*(SAB-Mf2(i-1,j)/MW2*SA*CB/SB2+Mf2(i,j)/MW2*(-CA)*SB/CB2)
	KLR(2,5,2,2,i,j) = hf_(i-1,j)*hf_(i,j)/sqrt2*CBA
	KLR(5,2,1,1,i,j) = KLR(2,5,1,1,i,j)
	KLR(5,2,2,2,i,j) = KLR(2,5,2,2,i,j)

	KLR(3,5,1,1,i,j) = cI*EL**2/(2d0*sqrt2*SW2)*(C2B-Mf2(i-1,j)/(MW2*TB2)+Mf2(i,j)/MW2*TB2) 
	KLR(3,5,2,2,i,j) = (0d0,0d0)
	KLR(5,3,1,1,i,j) = KLR(3,5,1,1,i,j)
	KLR(5,3,2,2,i,j) = KLR(3,5,2,2,i,j)

	KLR(4,5,1,1,i,j) = cI*EL**2/(2d0*sqrt2*SW2)*(S2B-Mf2(i-1,j)/(MW2*TB)-Mf2(i,j)/MW2*TB) 
	KLR(4,5,2,2,i,j) = cI*hf_(i-1,j)*hf_(i,j)/sqrt2
	KLR(5,4,1,1,i,j) = KLR(4,5,1,1,i,j)
	KLR(5,4,2,2,i,j) = KLR(4,5,2,2,i,j)

!neutral Higgs,Gp (Gm)

!error in the HiggsHunter p.408, necessary replacements: SA->-CA,CA->SA (wrong->correct)
	KLR(1,6,1,1,i,j) = EL**2/(2d0*sqrt2*SW2)*(-SAB-Mf2(i-1,j)/MW2*(-CA)/SB+Mf2(i,j)/MW2*SA/CB) 
	KLR(1,6,2,2,i,j) = -hf_(i-1,j)*hf_(i,j)/sqrt2*CBA
	KLR(6,1,1,1,i,j) = KLR(1,6,1,1,i,j)
	KLR(6,1,2,2,i,j) = KLR(1,6,2,2,i,j)

!error in the HiggsHunter p.408, necessary replacements: SA->CA,CA->-SA (wrong->correct)
	KLR(2,6,1,1,i,j) = EL**2/(2d0*sqrt2*SW2)*(CAB-Mf2(i-1,j)/MW2*(-SA)/SB-Mf2(i,j)/MW2*CA/CB)
	KLR(2,6,2,2,i,j) = hf_(i-1,j)*hf_(i,j)/sqrt2*SBA
	KLR(6,2,1,1,i,j) = KLR(2,6,1,1,i,j)
	KLR(6,2,2,2,i,j) = KLR(2,6,2,2,i,j)

	KLR(3,6,1,1,i,j) = cI*EL**2/(2d0*sqrt2*SW2)*(S2B-Mf2(i-1,j)/(MW2*TB)-Mf2(i,j)/MW2*TB) 
	KLR(3,6,2,2,i,j) =-cI*hf_(i-1,j)*hf_(i,j)/sqrt2
	KLR(6,3,1,1,i,j) = KLR(3,6,1,1,i,j)
	KLR(6,3,2,2,i,j) = KLR(3,6,2,2,i,j)

	KLR(4,6,1,1,i,j) =-cI*EL**2/(2d0*sqrt2*SW2)*(C2B+Mf2(i-1,j)/MW2-Mf2(i,j)/MW2) 
	KLR(4,6,2,2,i,j) = (0d0,0d0)
	KLR(6,4,1,1,i,j) = KLR(4,6,1,1,i,j)
	KLR(6,4,2,2,i,j) = KLR(4,6,2,2,i,j)


	endif
        
           enddo
        enddo

! neutral Higgs
	do a=1,6
        do b=1,6
	do i=1,2
	 do j=1,2
         do m=1,4
          do n=1,4
          do o=1,3

	  HHSfSf(a,b,i,j,m,n,o) = (0d0,0d0) 

	  do k=1,2
	   do l=1,2

            if(((a.lt.5).and.(b.lt.5)).or.((a.gt.4).and.(b.gt.4))) then	

	    HHSfSf(a,b,i,j,m,n,o) = HHSfSf(a,b,i,j,m,n,o) + KLR(a,b,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,o)*Krondelta(m,n)      !change coupling von L,R- to 1,2-basis

            elseif(((a.lt.5).and.(b.gt.4)).or.((a.gt.4).and.(b.lt.5))) then

              if((m.eq.1).or.(m.eq.3)) then

              HHSfSf(a,b,i,j,m,n,o) = HHSfSf(a,b,i,j,m,n,o) + KLR(a,b,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,o)*Krondelta(m+1,n)  !change coupling von L,R- to 1,2-basis
   
              else

              HHSfSf(a,b,i,j,m,n,o) = HHSfSf(a,b,i,j,m,n,o) + KLR(a,b,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,o)*Krondelta(m-1,n)  !change coupling von L,R- to 1,2-basis
 
              endif
 
            endif

   	   enddo
	  enddo
   
          enddo
          enddo
 	 enddo
         enddo
	enddo
	enddo
        enddo
        
!********************************************
! sFermion-sFermion-Vector-Vector-couplings *
!********************************************
!e.g. VVSfSf(xtype,ytype,sfmeigenup,sfmeigendown,sftypeup,sftypedown,sfgenup) no flavor mixing => sfgen2-independent, couplings multiplied by -i

!neglect flavor-mixing
!used Feynarts-coupling-files (MSSM- and MSSMQCD-modelfiles + Sparticles) 

!1=Gluon (=g), 2=Gamma (=G), 3=Z, 4=W+ (W-)


        do axhtype= 1,4                             !no flavor-mixing included, hence sfgen2 is not needed
          do ayhtype= 1,4
            do asfmeigen1= 1,2
              do asfmeigen2= 1,2
                do asftype1= 1,4
                  do asftype2=1,4
                    do asfgen1= 1,3

                    VVSfSf(axhtype,ayhtype,asfmeigen1,asfmeigen2,asftype1,asftype2,asfgen1)= (0d0,0d0)

                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do i= 1, 4                    
          do j=1, 4
            do k=1, 3
              
        do asfmeigen1=1,2
          do asfmeigen2=1,2

          if((i.eq.3).or.(i.eq.4)) then  

!g,g
       
           VVSfSf(1,1,asfmeigen1,asfmeigen2,i,j,k) = GS**2*Krondelta(asfmeigen1,asfmeigen2)*Krondelta(i,j)                                !see Sparticles p. 89; FeynArts MSSMQCD.dvi p.3            
  
           VVSfSft(1,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,1,asfmeigen1,asfmeigen2,i,j,k)
           VVSfSfu(1,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,1,asfmeigen1,asfmeigen2,i,j,k)
!G,g

	   VVSfSf(1,2,asfmeigen1,asfmeigen2,i,j,k) = 2d0*EL*GS*ef_(i)*Krondelta(asfmeigen1,asfmeigen2)*Krondelta(i,j)                     !see  FeynArts MSSMQCD.dvi p.3 + NeuQ2qx_couplings 
           VVSfSf(2,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,2,asfmeigen1,asfmeigen2,i,j,k)   
   
           VVSfSft(1,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,2,asfmeigen1,asfmeigen2,i,j,k)
           VVSfSfu(1,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,2,asfmeigen1,asfmeigen2,i,j,k)    
           VVSfSft(2,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,1,asfmeigen1,asfmeigen2,i,j,k)
           VVSfSfu(2,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,1,asfmeigen1,asfmeigen2,i,j,k)                                  

!g,Z

	   VVSfSf(1,3,asfmeigen1,asfmeigen2,i,j,k) = 2d0*EL*GS/(SW*CW)*(I3f_(i)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)-ef_(i)*SW2*Krondelta(asfmeigen1,asfmeigen2) )*Krondelta(i,j)  !see FeynArts MSSMQCD.dvi p.3
	   VVSfSf(3,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,3,asfmeigen1,asfmeigen2,i,j,k)   

           VVSfSft(1,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,3,asfmeigen1,asfmeigen2,i,j,k)
           VVSfSfu(1,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,3,asfmeigen1,asfmeigen2,i,j,k)    
           VVSfSft(3,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,1,asfmeigen1,asfmeigen2,i,j,k)
           VVSfSfu(3,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,1,asfmeigen1,asfmeigen2,i,j,k) 

!g,W+
           if(i.eq.3) then

            VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k) = sqrt2*EL*GS/SW*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i+1,j)        !see FeynArts MSSMQCD.dvi p.4
            VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k) 

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VVSfSft(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(1,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)  
            VVSfSft(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(4,1,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)      
            else
            VVSfSft(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k)
            endif
 
           else

            VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k) = sqrt2*EL*GS/SW*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i-1,j)        !see FeynArts MSSMQCD.dvi p.4
            VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k) 

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VVSfSft(1,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,1,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k)
            else
            VVSfSft(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(1,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(1,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(4,1,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,1,asfmeigen1,asfmeigen2,i,j,k)
            endif     

           endif  
 
          endif      

!GG

          VVSfSf(2,2,asfmeigen1,asfmeigen2,i,j,k) = 2d0*EL**2*ef_(i)**2*Krondelta(asfmeigen1,asfmeigen2)*Krondelta(i,j)                   !see Sparticles p. 512; EL=g2*SW, FeynArts MSSM.dvi p. 38

          VVSfSft(2,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,2,asfmeigen1,asfmeigen2,i,j,k)
          VVSfSfu(2,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,2,asfmeigen1,asfmeigen2,i,j,k)    

!G,Z
 
          VVSfSf(2,3,asfmeigen1,asfmeigen2,i,j,k) = 2d0*EL**2/(SW*CW)*(ef_(i)*I3f_(i)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)-ef_(i)**2*SW2*Krondelta(asfmeigen1,asfmeigen2) )*Krondelta(i,j) !see FeynArts MSSM.dvi p. 38 
	  VVSfSf(3,2,asfmeigen1,asfmeigen2,i,j,k) =  VVSfSf(2,3,asfmeigen1,asfmeigen2,i,j,k)      

          VVSfSft(2,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,3,asfmeigen1,asfmeigen2,i,j,k)
          VVSfSfu(2,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,3,asfmeigen1,asfmeigen2,i,j,k)    
          VVSfSft(3,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,2,asfmeigen1,asfmeigen2,i,j,k)
          VVSfSfu(3,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,2,asfmeigen1,asfmeigen2,i,j,k) 

!Z,Z
 
          VVSfSf(3,3,asfmeigen1,asfmeigen2,i,j,k) = 2d0*EL**2/(S2W**2)*((1d0 - 8*ef_(i)*I3f_(i)*SW2)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k) + 4*ef_(i)**2*SW2**2*Krondelta(asfmeigen1,asfmeigen2) )*Krondelta(i,j) !see Feynarts MSSM.dvi p. 39 + ???Sparticles p. 513???      

          VVSfSft(3,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,3,asfmeigen1,asfmeigen2,i,j,k)
          VVSfSfu(3,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,3,asfmeigen1,asfmeigen2,i,j,k)    

!W+,W-

          VVSfSf(4,4,asfmeigen1,asfmeigen2,i,j,k) = EL**2/(2*SW2)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i,j)             !see Feynarts MSSM.dvi p. 39 + Sparticles p. 513    

          VVSfSft(4,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,4,asfmeigen1,asfmeigen2,i,j,k)
          VVSfSfu(4,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,4,asfmeigen1,asfmeigen2,i,j,k)    

!W+-,G/Z

          if((i.eq.1).or.(i.eq.3)) then

            VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k) = EL**2/(3d0*sqrt2*SW)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i+1,j)  !see FeynArts MSSM.dvi p.39
            VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k) 

            VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k) = -EL**2/(3d0*sqrt2*CW)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i+1,j) !see FeynArts MSSM.dvi p.39
            VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k) 

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VVSfSft(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(2,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)  
            VVSfSft(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(4,2,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)   

            VVSfSft(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(3,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)  
            VVSfSft(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k) 
            VVSfSfu(4,3,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)         
            else
            VVSfSft(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k)

            VVSfSft(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k)   
            VVSfSfu(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k)
            endif


 
          else

            VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k) = EL**2/(3d0*sqrt2*SW)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i-1,j)  !see FeynArts MSSM.dvi p.39
            VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k) 

            VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k) = -EL**2/(3d0*sqrt2*CW)*USf(asfmeigen1,1,i,k)*USf(asfmeigen2,1,j,k)*Krondelta(i-1,j) !see FeynArts MSSM.dvi p.39
            VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k) 

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VVSfSft(2,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,2,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k)

            VVSfSft(3,4,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,3,asfmeigen1,asfmeigen2,i,j,k) = (0d0,0d0)    
            VVSfSfu(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k)
            else
            VVSfSft(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(2,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(2,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(4,2,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,2,asfmeigen1,asfmeigen2,i,j,k)

            VVSfSft(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(3,4,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(3,4,asfmeigen1,asfmeigen2,i,j,k)
            VVSfSft(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k)    
            VVSfSfu(4,3,asfmeigen1,asfmeigen2,i,j,k) = VVSfSf(4,3,asfmeigen1,asfmeigen2,i,j,k)
            endif     

          endif  


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
     &                    sqrt2*EL/SW*((ef_(b) - I3f_(b))*SW/CW*ZNeu(j,1) +  I3f_(b)*ZNeu(j,2))*USf(a,1,b,i))*Krondelta(b,m)*Krondelta(i,n)
     
                          AL(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*SB)*ZNeuC(j,4)*USf(a,1,b,i) + 
     &                    sqrt2*EL*ef_(b)/CW*ZNeuC(j,1)*USf(a,2,b,i))*Krondelta(b,m)*Krondelta(i,n)	

! --- incoming down-type sfermion---     

      		       else
      		       
                          AR(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*CB)*ZNeu(j,3)*USf(a,2,b,i) -
     &                    sqrt2*EL/SW*(I3f_(b)*ZNeu(j,2) + (ef_(b) - I3f_(b))*SW/CW*ZNeu(j,1))*USf(a,1,b,i))*Krondelta(b,m)*Krondelta(i,n)

                          AL(j,a,b,m,i,n) =  (- EL*Mf(b,i)/(SW*sqrt2*MW*CB)*ZNeuC(j,3)*USf(a,1,b,i) + 
     &                    sqrt2*EL*ef_(b)/CW*ZNeuC(j,1)*USf(a,2,b,i))*Krondelta(b,m)*Krondelta(i,n)
                                 
	               endif
	          enddo

	          
! ------------------------------ Vertex Chargino - Sfermion (in) - Fermion (out) ------------------------------                                                         !differs from couplings implemented in CalcHEP->not comparable

	          do j= 5, 6   

! --- incoming up-type sfermion ---	          

	          if ((b.eq.1).or.(b.eq.3)) then

! 	         	  KL(j,a,b,i) =  -EL/SW*dconjg(UCha(j-4,1))*dconjg(USf(a,1,b+1,i)) + EL*Mf(b+1,i)/(SW*sqrt2*MW*CB)*dconjg(UCha(j-4,2))*dconjg(USf(a,2,b+1,i))   != C_L in Sparticles p. 515
! 	         	  KR(j,a,b,i) =  EL*Mf(b,i)/(SW*sqrt2*MW*SB)*VCha(j-4,2)*dconjg(USf(a,1,b+1,i))                                                                 != F_R in Sparticles p. 515 
 	         	
	         	  AL(j,a,b,m,i,n) =  EL*Mf(b+1,i)/(SW*sqrt2*MW*CB)*dconjg(UCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m-1)*Krondelta(i,n)                             != E_R in Sparticles p. 515 , no t_i-e-coupl. etc. and no generation-mixing                                                                                                                                                                        
 	         	  AR(j,a,b,m,i,n) = (-EL/SW*VCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)/(SW*sqrt2*MW*SB)*VCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m-1)*Krondelta(i,n)    != D_L in Sparticles p. 515

! --- incoming down-type sfermion --- 
                        
 	         	else

! 	         	  KL(j,a,b,i) =  -EL/SW*dconjg(VCha(j-4,1))*dconjg(USf(a,1,b-1,i)) + EL*Mf(b-1,i)/(SW*sqrt2*MW*SB)*dconjg(VCha(j-4,2))*dconjg(USf(a,2,b-1,i))   != D_L in Sparticles p. 515
! 	         	  KR(j,a,b,i) =  EL*Mf(b,i)/(SW*sqrt2*MW*CB)*UCha(j-4,2)*dconjg(USf(a,1,b-1,i))                                                                 != E_R in Sparticles p. 515

 	         	  AL(j,a,b,m,i,n) = EL*Mf(b-1,i)/(SW*sqrt2*MW*SB)*dconjg(VCha(j-4,2))*USf(a,1,b,i)*Krondelta(b,m+1)*Krondelta(i,n)                              != F_R in Sparticles p. 515
 	         	  AR(j,a,b,m,i,n) = (-EL/SW*UCha(j-4,1)*USf(a,1,b,i) + EL*Mf(b,i)/(SW*sqrt2*MW*CB)*UCha(j-4,2)*USf(a,2,b,i))*Krondelta(b,m+1)*Krondelta(i,n)    != C_L in Sparticles p. 515 
 	         	
 	         	endif
	          enddo

                  enddo
                enddo
 	      enddo
	   enddo
        enddo
	 
!******************************************************** 
! Vertex Gluino  - Squark (outgoing) - Quark (incoming) *
!********************************************************
!  e.g. GluL/R(sfmeigenin,sftypein,ftypeout,sfgenin,fgenout), couplings multiplied by -i

           do a=1, 2                                                                              !taken from NeuQ2qx + Sparticles p. 88 
             do b= 1, 4
               do m= 1, 4
                 do i= 1, 3
                    do n= 1, 3

                     if((b.gt.2).and.(m.gt.2)) then

                        GluL(a,b,m,i,n) = - GS*sqrt2*USf(a,1,b,i)*Krondelta(b,m)*Krondelta(i,n)
                        GluR(a,b,m,i,n) =   GS*sqrt2*USf(a,2,b,i)*Krondelta(b,m)*Krondelta(i,n)

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

c	C(1,3) = -findValW("SC")
c	C(1,4) = -findValW("SC")

	C(1,3) = -GS
	C(1,4) = -GS

	C(2,1) = -ef_(1)*EL
	C(2,2) = -ef_(2)*EL
	C(2,3) = -ef_(3)*EL
	C(2,4) = -ef_(4)*EL

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
     
     
           VSfSf(3,i,j,k,l,m,n)=-EL/(SW*CW)*(I3f_(k)*USf(j,1,l,n)*USf(i,1,k,m)-SW2*ef_(k)*Krondelta(j,i))*Krondelta(k,l)*Krondelta(m,n)  !Z-Boson
           
           VSfSft(3,i,j,k,l,m,n) = VSfSf(3,i,j,k,l,m,n) 
           VSfSfu(3,i,j,k,l,m,n) = VSfSf(3,i,j,k,l,m,n)           
     
! --- incoming up-type sfermion ---
               
           if ((l.eq.1).or.(l.eq.3)) then                                                                                                !W+/- -Boson

    	    VSfSf(4,i,j,k,l,m,n) = - EL/(sqrt2*SW)*USf(j,1,l,n)*USf(i,1,k,m)*Krondelta(l+1,k)*Krondelta(m,n)                          

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VSfSft(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n) 
            VSfSfu(4,i,j,k,l,m,n) = (0d0,0d0)    
            else
            VSfSft(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n)   
            VSfSfu(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n)
            endif

! --- incoming down-type sfermion ---
               
           else
    	       
    	    VSfSf(4,i,j,k,l,m,n) = - EL/(sqrt2*SW)*USf(j,1,l,n)*USf(i,1,k,m)*Krondelta(l-1,k)*Krondelta(m,n)         

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then                                                   !u-channel only needs to be switched of for '+-' final-states (e.g. W+H- etc.)
            VSfSft(4,i,j,k,l,m,n) = (0d0,0d0)    
            VSfSfu(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n)
            else
            VSfSft(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n)    
            VSfSfu(4,i,j,k,l,m,n) = VSfSf(4,i,j,k,l,m,n)
            endif
                                
           endif
               

                    enddo
                  enddo 
                enddo
 	      enddo
 	    enddo
	  enddo

!**************************
! Higgs-Sfermion-Sfermion *
!**************************
!e.g. HSfSf(xhtype,sfmeigenin,sfmeigenout,sftypein,sftypeout,sfgenin,sfgenout), couplings multiplied by -i    

!From Higgs-Hunter (p.394 f)

        do i= 1, 4
           do j=1, 3


	if ((i.eq.1).or.(i.eq.3)) then

	GLR(1,1,1,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*CA + EL*MZ*(I3f_(i)-ef_(i)*SW2)*SAB/(SW*CW)             
	GLR(1,2,2,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*CA + EL*MZ*ef_(i)*SW2*SAB/(SW*CW)                      
	GLR(1,2,1,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*CA + MUE*SA) 
	GLR(1,1,2,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*CA + MUE*SA) 

	GLR(2,1,1,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*SA - EL*MZ*(I3f_(i)-ef_(i)*SW2)*CAB/(SW*CW)             
	GLR(2,2,2,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*SA - EL*MZ*ef_(i)*SW2*CAB/(SW*CW)                       
	GLR(2,2,1,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*SA - MUE*CA)                                           
	GLR(2,1,2,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*SA - MUE*CA)                                           
      
	GLR(3,1,1,i,j) = (0d0,0d0)
	GLR(3,2,2,i,j) = (0d0,0d0)
	GLR(3,2,1,i,j) = cI*hf_(i,j)/sqrt2*(Af(i,j)*CB + MUE*SB) 
	GLR(3,1,2,i,j) = -cI*hf_(i,j)/sqrt2*(Af(i,j)*CB + MUE*SB) 

	GLR(4,1,1,i,j) = (0d0,0d0)
	GLR(4,2,2,i,j) = (0d0,0d0)
	GLR(4,2,1,i,j) = cI*hf_(i,j)/sqrt2*(Af(i,j)*SB - MUE*CB) 
	GLR(4,1,2,i,j) = -cI*hf_(i,j)/sqrt2*(Af(i,j)*SB - MUE*CB) 

! G_{LR}^{tb} where i = up-type
	GLR(5,1,1,i,j) = hf_(i+1,j)*Mf(i+1,j)*SB + hf_(i,j)*Mf(i,j)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,i,j) = hf_(i,j)*Mf(i+1,j)*CB + hf_(i+1,j)*Mf(i,j)*SB 
	GLR(5,2,1,i,j) = hf_(i,j)*(Af(i,j)*CB + MUE*SB)
	GLR(5,1,2,i,j) = hf_(i+1,j)*(Af(i+1,j)*SB + MUE*CB)

	GLR(6,1,1,i,j) = -hf_(i+1,j)*Mf(i+1,j)*CB + hf_(i,j)*Mf(i,j)*SB + EL*MW/(sqrt2*SW)*C2B
	GLR(6,2,2,i,j) = hf_(i,j)*Mf(i+1,j)*SB - hf_(i+1,j)*Mf(i,j)*CB                                  
	GLR(6,2,1,i,j) = hf_(i,j)*(Af(i,j)*SB - MUE*CB)
	GLR(6,1,2,i,j) = hf_(i+1,j)*(-Af(i+1,j)*CB + MUE*SB)

	else

	GLR(1,1,1,i,j) = sqrt2*hf_(i,j)*Mf(i,j)*SA + EL*MZ*(I3f_(i)-ef_(i)*SW2)*SAB/(SW*CW) 
	GLR(1,2,2,i,j) = sqrt2*hf_(i,j)*Mf(i,j)*SA + EL*MZ*ef_(i)*SW2*SAB/(SW*CW) 
	GLR(1,2,1,i,j) = hf_(i,j)/sqrt2*(Af(i,j)*SA + MUE*CA) 
	GLR(1,1,2,i,j) = hf_(i,j)/sqrt2*(Af(i,j)*SA + MUE*CA) 

	GLR(2,1,1,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*CA - EL*MZ*(I3f_(i)-ef_(i)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,i,j) = -sqrt2*hf_(i,j)*Mf(i,j)*CA - EL*MZ*ef_(i)*SW2*CAB/(SW*CW) 
	GLR(2,2,1,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*CA - MUE*SA) 
	GLR(2,1,2,i,j) = -hf_(i,j)/sqrt2*(Af(i,j)*CA - MUE*SA) 

	GLR(3,1,1,i,j) = (0d0,0d0)
	GLR(3,2,2,i,j) = (0d0,0d0)
	GLR(3,2,1,i,j) = cI*hf_(i,j)/sqrt2*(Af(i,j)*SB + MUE*CB) 
	GLR(3,1,2,i,j) = -cI*hf_(i,j)/sqrt2*(Af(i,j)*SB + MUE*CB) 

	GLR(4,1,1,i,j) = (0d0,0d0)
	GLR(4,2,2,i,j) = (0d0,0d0)
	GLR(4,2,1,i,j) = -cI*hf_(i,j)/sqrt2*(Af(i,j)*CB - MUE*SB) 
	GLR(4,1,2,i,j) = cI*hf_(i,j)/sqrt2*(Af(i,j)*CB - MUE*SB) 

! G_{LR}^{bt} where i = down-type
	GLR(5,1,1,i,j) = hf_(i,j)*Mf(i,j)*SB + hf_(i-1,j)*Mf(i-1,j)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,i,j) = hf_(i-1,j)*Mf(i,j)*CB + hf_(i,j)*Mf(i-1,j)*SB
	GLR(5,2,1,i,j) = hf_(i,j)*(Af(i,j)*SB + MUE*CB)
        GLR(5,1,2,i,j) = hf_(i-1,j)*(Af(i-1,j)*CB + MUE*SB) 

	GLR(6,1,1,i,j) = -hf_(i,j)*Mf(i,j)*CB + hf_(i-1,j)*Mf(i-1,j)*SB + EL*MW/(sqrt2*SW)*C2B
	GLR(6,2,2,i,j) = hf_(i-1,j)*Mf(i,j)*SB - hf_(i,j)*Mf(i-1,j)*CB                                 
	GLR(6,2,1,i,j) = hf_(i,j)*(-Af(i,j)*CB + MUE*SB)
	GLR(6,1,2,i,j) = hf_(i-1,j)*(Af(i-1,j)*SB - MUE*CB) 

	endif
        
           enddo
        enddo

! neutral Higgs

	do a=1,4
	do i=1,2
	 do j=1,2
         do m=1,4
          do n=1,4
          do o=1,3
            do p=1,3

	  HSfSf(a,i,j,m,n,o,p) = (0d0,0d0) 

	  do k=1,2
	   do l=1,2
	
	    HSfSf(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p) + GLR(a,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,p)*Krondelta(m,n)*Krondelta(o,p)

   	   enddo
	  enddo
   
           
            HSfSft(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)                
            HSfSfu(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
    
           enddo
          enddo
          enddo
 	 enddo
         enddo
	enddo
	enddo

! charged Higgs

	do a=5,6
	do i=1,2
	 do j=1,2
         do m=1,4                                 
          do n=1,4
          do o=1,3
           do p=1,3
	
	  HSfSf(a,i,j,m,n,o,p) = (0d0,0d0) 
	
	  do k=1,2
	   do l=1,2
	
         if((m.eq.1).or.(m.eq.3))  then   !up-case

	    HSfSf(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p) + GLR(a,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,p)*Krondelta(m+1,n)*Krondelta(o,p)
         
         else                             !down-case

            HSfSf(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p) + GLR(a,k,l,m,o)*USf(i,k,m,o)*USf(j,l,n,p)*Krondelta(m-1,n)*Krondelta(o,p) 

         endif

	   enddo
	  enddo          

         if((m.eq.1).or.(m.eq.3))  then 

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then  !u-channel only needs to be switched of for +- final-states (e.g. W+H- etc.)
            HSfSft(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)                                
            HSfSfu(a,i,j,m,n,o,p) = (0d0,0d0)
            else
            HSfSft(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
            HSfSfu(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
            endif
 
         else

            if(((xtype.eq.4).or.(xtype.eq.8)).and.((ytype.eq.4).or.(ytype.eq.8))) then  !u-channel only needs to be switched of for +- final-states (e.g. W+H- etc.)
            HSfSft(a,i,j,m,n,o,p) = (0d0,0d0)
            HSfSfu(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
            else
            HSfSft(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
            HSfSfu(a,i,j,m,n,o,p) = HSfSf(a,i,j,m,n,o,p)
            endif

         endif

           enddo
          enddo
          enddo
         enddo
	 enddo
	enddo
	enddo

! *****************************************************************************************************               
! *                               Vertex Squark - Squark - Squark - Squark                            *
! *****************************************************************************************************
	


!aus Julias code: derzeitig nur Anteile der 4-sQuark-Kopplung die prop. zu GS sind implementiert, Vergleiche auch die hier implementierte Kopplung mit z.B. Sparticles S. 215ff
!CrossCheck mit meinem sF4a und SF4b ergab Übereinstimmung (*4 zu Julias Kopplung für Farbe)

!	do i=1,2
!	   do j=1,2
!        	do k=1,2
!              	   do l=1,2
!		       do m = 1, 4
!		          do n = 1, 4
!		             do o = 1, 3

!              sF4(i,k,j,l,m,n,o) = -GS**2*(USf(j,1,m,o)*USf(k,1,m,o)-USf(j,2,m,o)*USf(k,2,m,o))*
!     &              (USf(l,1,n,o)*USf(i,1,n,o)-USf(l,2,n,o)*USf(i,2,n,o))


!			      enddo
!			  enddo
!		      enddo
!		  enddo
!	      enddo
!	   enddo
!     	enddo

!from Rosiek p. 47 crosschecked with Sparticles p.520  (only strong-coupling-part (~GS^2) implemented yet, no generation-mixing, couplings multiplied by -i)
!until now only checked for squarks 

!These Couplings reduce to Julias coupling (*colorfactor) above in the case of sQuark-Loops
!Convention i,k=>in , l,j=> out 
        do i = 1,2
	   do j = 1,2
              do k = 1,2
              	  do l = 1,2
		      do m = 1,4
		          do n = 1,4
		              do o = 1,3
      
!general structure colora*sf4a+colorb*sF4b with colorfactors from Rosiek p. 48
        
              if(m.eq.n) then

!e.g. sF4(sfmeigenin1,sfmeigenin2,sfmeigenout1,sfmeigenout2,sftype1in,sftype2in,sfgenin)

              Xil = Krondelta(i,l)- 2d0*USf(i,1,m,o)*USf(l,1,n,o)
              Xkj = Krondelta(k,j)- 2d0*USf(k,1,n,o)*USf(j,1,m,o)
              Xij = Krondelta(i,j)- 2d0*USf(i,1,m,o)*USf(j,1,m,o)
              Xkl = Krondelta(k,l)- 2d0*USf(k,1,n,o)*USf(l,1,n,o)

              sF4a(i,k,j,l,m,n,o) = -GS**2/6d0*(3*Xil*Xkj-Xij*Xkl)    
              sF4b(i,k,j,l,m,n,o) = -GS**2/6d0*(3*Xij*Xkl-Xil*Xkj)

              else

!e.g. sF4(sfmeigenin1up,sfmeigenin1down,sfmeigenout2up,sfmeigenout2down,sftypeinup,sftypeindown,sfgenin)

              Xij = Krondelta(i,j)-2d0*USf(i,1,m,o)*USf(j,1,m,o)
              Xkl = Krondelta(k,l)-2d0*USf(k,1,n,o)*USf(l,1,n,o)

              sF4a(i,k,j,l,m,n,o) = -GS**2/2d0*Xij*Xkl                 
              sF4b(i,k,j,l,m,n,o) =  GS**2/6d0*Xij*Xkl

              endif

			      enddo
			  enddo
		      enddo
		  enddo
	      enddo
	   enddo
     	enddo


!****************************************************************************************************************

	end


!functions used in couplings

! ===============================================================================================================  	
! ============================================== Fermion charges ================================================
! ===============================================================================================================


	double precision function ef_(type)
	integer type
	
	if (type.eq.1) then
	
	 ef_ = 0d0
	 
	elseif (type.eq.2) then
	
	 ef_ = -1d0
	 
	elseif (type.eq.3) then
	
	 ef_ = 2/3d0
	 
	elseif (type.eq.4) then
	
	 ef_ = -1/3d0
	 
	else 
	endif
	end



	
! ===============================================================================================================
! ============================================== Fermion isospin ================================================
! ===============================================================================================================

	double precision function I3f_(type)
	integer type
	
	if ((type.eq.1).or.(type.eq.3)) then
	
  	 I3f_ = 1/2d0
  	 
	 else
	 
  	 I3f_ = -1/2d0
  	 
	endif	
	end



	
! ===============================================================================================================	
! ============================================== Yukawa coupling ================================================
! ===============================================================================================================

	function hf_(type, gen)
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Model.h"
	integer type, gen
	
	if ((type.eq.1).or.(type.eq.3)) then
	
  	 hf_ = EL*Mf(type,gen)/(sqrt2*MW*SW*SB) 
  	 
	else
	
  	 hf_ = EL*Mf(type,gen)/(sqrt2*MW*SW*CB) 
  	 
	endif	
	end
	
	
