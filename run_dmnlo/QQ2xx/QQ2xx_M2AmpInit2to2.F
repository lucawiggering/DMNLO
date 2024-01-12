************************************************

	subroutine QQ2xx_M2AmpInit2to2(mytree,channel_,corrProp_,corrVert_,corrBox_,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)

	implicit none

#include "QQ2xx_Model.h"

         integer mytree,channel_,corrProp_,corrVert_,corrBox_
         double precision QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box,AmpTree,AmpProp,AmpVert,AmpBox,Weight

         external QQ2xx_Init2to2,QQ2xx_M2AmpHH,QQ2xx_M2AmpVH QQ2xx_M2AmpHV,QQ2xx_M2AmpVV,QQ2xx_M2AmpFFbar,QQ2xx_M2AmpFF,QQ2xx_M2onshellsubtract

! initialize flags and switches for propagators, colorfactors etc.

         call QQ2xx_Init2to2(xtype,ytype)                                                                                                         

! initialize O(alpha_s) corrections (in QQ2xx_M2AmpCorr.F)

         call QQ2xx_CorrSX(mytree)
         call QQ2xx_CorrVX(mytree)
         call QQ2xx_CorrTXProp(mytree)
         call QQ2xx_CorrUXProp(mytree)

         QQ2xx_M2Tree = 0d0
         QQ2xx_M2Prop = 0d0
         QQ2xx_M2Vert = 0d0
         QQ2xx_M2Box = 0d0

! All vector boson final states get the Goldstone contributions added i.e. ZZ + ZG + GZ + GG
! Higgs boson final states untouched
! Fermion final states qq vs qqbar differentiation
	 if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.lt.5).and.(ytype.lt.5))) then

                      call QQ2xx_M2onshellsubtract(xtype,ytype,0,0,Weight)
                      call QQ2xx_M2AmpVV(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,ytype,AmpTree,AmpProp,AmpVert,AmpBox)
	 	
                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                                                
	 	
                if((xtype.eq.3).and.(ytype.eq.3)) then                                                             !add Goldstones for Z and W+/- final-state 
    
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,4,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox
                                      

                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,4,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox
                   

                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,4,4,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                

                elseif((xtype.eq.4).and.(ytype.eq.4)) then    
   
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,6,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox


                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,6,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox

 
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,6,6,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  
 
                elseif((xtype.eq.4).and.(ytype.eq.3)) then         

                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,4,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox     
                                         

                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,6,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox

 
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,6,4,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  

                elseif((xtype.eq.3).and.(ytype.eq.4)) then         
 
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,6,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                                             


                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,4,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox

 
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,4,6,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  

                elseif((xtype.lt.3).and.(ytype.eq.3)) then                                                         !for GZ/AZ final-state

                      call QQ2xx_M2onshellsubtract(xtype,0,0,4,Weight)
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,4,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox   
                                      
                elseif((xtype.eq.3).and.(ytype.lt.3)) then                                                         !for ZG/ZA final state

                      call QQ2xx_M2onshellsubtract(0,ytype,4,0,Weight)
                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,4,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox          
                
                elseif((xtype.lt.3).and.(ytype.eq.4)) then                                                         !for GW-/AW- final-state

                      call QQ2xx_M2onshellsubtract(xtype,0,0,6,Weight)
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,6,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  
                                      
                elseif((xtype.eq.4).and.(ytype.lt.3)) then                                                         !for W+G/W+A final state

                      call QQ2xx_M2onshellsubtract(0,ytype,6,0,Weight)
                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,6,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                                                                                 
                                                                                                                   
                endif

!                print*,'VV'
	 	
	 else if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.gt.4).and.(ytype.gt.4))) then
	 	
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,xhtype,yhtype,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  
	 	
!                 print*,'HH'

         else if (((f3type.eq.0).and.(f4type.eq.0)).and.((ytype.gt.4).and.(xtype.lt.5))) then
		
                      call QQ2xx_M2onshellsubtract(xtype,0,0,yhtype,Weight)
                      call QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,xtype,yhtype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  
	 	
                if(xtype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 3 = vector, particle 4 = scalar)

                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,4,yhtype,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox         

                elseif((xtype.eq.4)) then   
                        
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,6,yhtype,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                                      

                endif

!                print*,'VH'

         else if (((f3type.eq.0).and.(f4type.eq.0)).and.((xtype.gt.4).and.(ytype.lt.5))) then			

                      call QQ2xx_M2onshellsubtract(0,ytype,xhtype,0,Weight)
                      call QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,xhtype,ytype,AmpTree,AmpProp,AmpVert,AmpBox)

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree - Weight*AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox

                if(ytype.eq.3) then                                                                                !add Goldstones for Z and W+ final-state (here particle 4 = vector, particle 3 = scalar)

                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,xhtype,4,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox       

                elseif((ytype.eq.4)) then                      
       
                      call QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,xhtype,6,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox                                 

                endif

!                print*,'HV'

	 else if ((ytype.eq.0).and.(xtype.eq.0).and.(f3type.ne.0).and.(f4type.ne.0)) then  

                if(n3*n4.lt.0) then                                                                                            

                      call QQ2xx_M2AmpFFbar(mytree,channel_,corrVert_,f3type,f3gen,f4type,f4gen,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox    

!                      print*,'FFbartree'
 
                else

                      call QQ2xx_M2AmpFF(channel_,f3type,f3gen,f4type,f4gen,AmpTree,AmpProp,AmpVert,AmpBox) 

                      QQ2xx_M2Tree = QQ2xx_M2Tree + AmpTree  

                      QQ2xx_M2Prop = QQ2xx_M2Prop + AmpProp

                      QQ2xx_M2Vert = QQ2xx_M2Vert + AmpVert

                      QQ2xx_M2Box = QQ2xx_M2Box + AmpBox  

!                      print*,'FFtree'

                endif

         else

	        QQ2xx_M2Tree = 0d0

                print*,'Final-state not covered by QQ2xx_Tree(2->2)'

	 endif

	
	end
