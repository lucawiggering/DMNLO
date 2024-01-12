************************************************

	 double precision function QQ2xx_M2Tree23()

	 implicit none

#include "QQ2xx_Model.h"	

	 double precision QQ2xx_M2Tree23VV, QQ2xx_M2Tree23VH, QQ2xx_M2Tree23HV, QQ2xx_M2Tree23HH, QQ2xx_M2Tree23FFbar, QQ2xx_M2Tree23FF

         external QQ2xx_Init2to3

!initialize propagators, colorfactors etc.

         call QQ2xx_Init2to3()                                                                                     

         QQ2xx_M2Tree23 = 0d0


	 if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.lt.5).and.(ytype.lt.5))) then
	 	
                QQ2xx_M2Tree23 = QQ2xx_M2Tree23VV()                                                               
	 	
                if((xtype.eq.3).and.(ytype.eq.3)) then                                                             !add Goldstones for Z and W+/- final-state 

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,4)                                        

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(4,ytype)                                        

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(4,4)       

                elseif((xtype.eq.4).and.(ytype.eq.4)) then    

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,6) 

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(6,ytype) 
 
                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(6,6)
 
                elseif((xtype.eq.4).and.(ytype.eq.3)) then                                                       

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,4) 

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(6,ytype) 
 
                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(6,4)
                
                elseif((xtype.eq.3).and.(ytype.eq.4)) then                                                       

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,6) 

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(4,ytype) 
 
                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(4,6)

                elseif((xtype.lt.3).and.(ytype.eq.3)) then                                                         !for GZ/AZ final-state

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,4)  
                                      
                elseif((xtype.eq.3).and.(ytype.lt.3)) then                                                         !for ZG/ZA final state

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(4,ytype)             
                
                elseif((xtype.lt.3).and.(ytype.eq.4)) then                                                         !for GW-/AW- final-state

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23VH(xtype,6)  
                                      
                elseif((xtype.eq.4).and.(ytype.lt.3)) then                                                         !for W+G/W+A final state

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HV(6,ytype)                                                                                             
                                                                                                                   
                endif

!                print*,'23VV'
	 	
	 else if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.gt.4).and.(ytype.gt.4))) then
	 	
	 	QQ2xx_M2Tree23 = QQ2xx_M2Tree23HH(xhtype,yhtype)
	 	
!                print*,'23HH'

	 else if (((f3type.eq.0).and.(f4type.eq.0)).and.((ytype.gt.4).and.(xtype.lt.5))) then
		
                QQ2xx_M2Tree23 = QQ2xx_M2Tree23VH(xtype,yhtype)
	 	
                if(xtype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 3 = vector, particle 4 = scalar)

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(4,yhtype)   

                elseif((xtype.eq.4)) then   

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(6,yhtype)                                 !G+H-

                endif

!                print*,'23VH'

         else if (((f3type.eq.0).and.(f4type.eq.0)).and.((ytype.lt.5).and.(xtype.gt.4))) then
		
                QQ2xx_M2Tree23 = QQ2xx_M2Tree23HV(xhtype,ytype)

                if(ytype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 4 = vector, particle 3 = scalar)

                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(xhtype,4)  

                elseif((ytype.eq.4)) then                      
       
                      QQ2xx_M2Tree23 = QQ2xx_M2Tree23 + QQ2xx_M2Tree23HH(xhtype,6)                                 !H+G- 

                endif

!                print*,'23HV'

	 else if ((ytype.eq.0).and.(xtype.eq.0).and.(f3type.ne.0).and.(f4type.ne.0)) then 

         	      QQ2xx_M2Tree23 = QQ2xx_M2Tree23FFbar()

!               print*,'FFbar'

         else

	        QQ2xx_M2Tree23 = 0d0

                print*,'Final-state not covered by QQ2xx_Tree(2->3)'

	 endif

	
	end
