************************************************

      subroutine QQ2xx_M2DigestAmp(costh)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"

      double precision costh,kappa

      external QQ2xx_M2AmpDigestVV, QQ2xx_M2AmpDigestVH, QQ2xx_M2AmpDigestHV, QQ2xx_M2AmpDigestHH, QQ2xx_M2AmpDigestFF, QQ2xx_M2AmpDigestFFbar

! Mandelstam variables t=(p1-k1)^2 , u=(p1-k2)^2
      tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &       0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &       0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

      uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
      ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

!initialize O(alpha_s) corrections
      call QQ2xx_CorrSX(0)
      call QQ2xx_CorrVX(0)
      call QQ2xx_CorrTXProp(0)
      call QQ2xx_CorrUXProp(0)         

      if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.lt.5).and.(ytype.lt.5))) then
	 	
                        call QQ2xx_M2AmpDigestVV(xtype,ytype)                                                                    
	 	
                if((xtype.eq.3).and.(ytype.eq.3)) then                                                             !add Goldstones for Z and W+/- final-state 
 
                        call QQ2xx_M2AmpDigestVH(xtype,4)                                         
  
                        call QQ2xx_M2AmpDigestHV(4,ytype)                                   
 
                        call QQ2xx_M2AmpDigestHH(4,4)  

                elseif((xtype.eq.4).and.(ytype.eq.4)) then    
  
                        call QQ2xx_M2AmpDigestVH(xtype,6)   

                        call QQ2xx_M2AmpDigestHV(6,ytype) 

                        call QQ2xx_M2AmpDigestHH(6,6) 
 
                elseif((xtype.eq.4).and.(ytype.eq.3)) then                                                       

                        call QQ2xx_M2AmpDigestVH(xtype,4)   

                        call QQ2xx_M2AmpDigestHV(6,ytype)  

                        call QQ2xx_M2AmpDigestHH(6,4) 

                elseif((xtype.eq.3).and.(ytype.eq.4)) then                                                       
  
                        call QQ2xx_M2AmpDigestVH(xtype,6)   

                        call QQ2xx_M2AmpDigestHV(4,ytype) 

                        call QQ2xx_M2AmpDigestHH(4,6) 

                elseif((xtype.lt.3).and.(ytype.eq.3)) then                                                         !for GZ/AZ final-state
  
                        call QQ2xx_M2AmpDigestVH(xtype,4)    
                                      
                elseif((xtype.eq.3).and.(ytype.lt.3)) then                                                         !for ZG/ZA final state

                        call QQ2xx_M2AmpDigestHV(4,ytype)              
                
                elseif((xtype.lt.3).and.(ytype.eq.4)) then                                                         !for GW-/AW- final-state
   
                        call QQ2xx_M2AmpDigestVH(xtype,6)    
                                      
                elseif((xtype.eq.4).and.(ytype.lt.3)) then                                                         !for W+G/W+A final state

                        call QQ2xx_M2AmpDigestHV(6,ytype)                                                                                          
                                                                                                                   
                endif
	 	
      else if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.gt.4).and.(ytype.gt.4))) then

                        call QQ2xx_M2AmpDigestHH(xhtype,yhtype) 

      else if (((f3type.eq.0).and.(f4type.eq.0)).and.((ytype.gt.4).and.(xtype.lt.5))) then

                        call QQ2xx_M2AmpDigestVH(xtype,yhtype)  
	 	
                if(xtype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 3 = vector, particle 4 = scalar)

                        call QQ2xx_M2AmpDigestHH(4,yhtype)   
         

                elseif((xtype.eq.4)) then   

                        call QQ2xx_M2AmpDigestHH(6,yhtype)  

                endif

!                print*,'VH'

      else if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.gt.4).and.(ytype.lt.5))) then			

                        call QQ2xx_M2AmpDigestHV(xhtype,ytype) 

                if(ytype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 4 = vector, particle 3 = scalar)

                        call QQ2xx_M2AmpDigestHH(xhtype,4)   

                elseif((ytype.eq.4)) then                      

                        call QQ2xx_M2AmpDigestHH(xhtype,6) 

                endif

!                print*,'HV'

      else if ((ytype.eq.0).and.(xtype.eq.0).and.(f3type.ne.0).and.(f4type.ne.0)) then  

                if(n3*n4.lt.0) then                                                                                !fermion-antifermion-case             

                        call QQ2xx_M2AmpDigestFFbar(f3type,f3gen,f4type,f4gen) 

!                      print*,'FFbartree'
 
                else

                        call QQ2xx_M2AmpDigestFF(f3type,f3gen,f4type,f4gen)

!                      print*,'FFtree'

                endif

      else

                print*,'Final-state not covered by QQ2xx_Tree(2->2)'

      endif

	
      end

