!**********************************************************************************************************************************


	subroutine QQ2xx_Init2to2(x_arg,y_arg)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_Init2to2.h"	

	integer i,x_arg,y_arg

!==================================================================================================================================
!                                                      ----Propagators----
!==================================================================================================================================

	DenSmV(1)  =  sqrtS**2                                                     !Gluon

	DenSmV(2)  =  sqrtS**2                                                     !Photon

	DenSmV(3)  =  (sqrtS**2 - MZ2) + cI*MZ*wZ                                  !Z-Boson

	DenSmV(4)  =  (sqrtS**2 - MW2) + cI*MW*wW                                  !W-Boson

!=============================----Index:1=h,2=H,3=P,4=G0,5=H+,6=G+----=============================================================
 
        do i = sumstartH,sumendH                                                            
                                                             
	  DenSmH(i) = (sqrtS**2 - MHk(i)**2) + cI*MHk(i)*wHk(i)                    !Higgs-Boson

	enddo

!================================================================================================================================
!                                          ----Gluon/Gluino-exchange (QCD-corrections)----
!================================================================================================================================

                                                       Looptype0 = 1                     
                                                       Loopmeigen1 = sfmeigen1                   
                                                       Looptype1 = sftype1
                                                       Loopgen1 = sfgen1
                                                       Loopmeigen2 = sfmeigen2
                                                       Looptype2 = sftype2
                                                       Loopgen2 = sfgen2

                                                       Looptype0 = 1                    
                                                       Loopmeigen1 = sfmeigen1
                                                       Looptype1 = sftype1
                                                       Loopgen1 = sfgen1
                                                       Loopmeigen2 = sfmeigen2
                                                       Looptype2 = sftype2
                                                       Loopgen2 = sfgen2

                                                       Looptypet = proptypet
                                                       Loopgent = propgent

                                                       Looptypeu = proptypeu
                                                       Loopgenu = propgenu

!needed in propagator corrections

                                                       Looptypeout1 = proptypet
                                                       Loopgenout1 = propgent      !No generation-mixing
                                                       Looptypeout2= proptypeu
                                                       Loopgenout2 = propgenu      !No generation-mixing

!==================================================================================================================================
!                                                      ----colorfactors----                                                         
!==================================================================================================================================
!Note that for strong final states at NLO, this implementation would lead to a double counting of colorfactors since 
!colorfactors would then be included in the tree level AND in the higher order corrections.

!colorfactors appearing in higher order corrections
                        if(Looptype0.eq.1) then                                    !SUSY-QCD corrections  (standard case)              
                                           
                                                       colorfactorSSVV = CF      
                                                       colorfactorSSVGl = CF 
                                                       colorfactorSSSV = CF 
                                                       colorfactorSSSSV = CF
                                                       colorfactorSSSGl = CF 
                                                       colorfactorQQVV = CF
                                                       colorfactorQQQQ = 1d0       !colorfactor already included in coupling

                                                       colorfactorVS = CF    
                                                       colorfactorV = 0d0          !in the case of Gluons taken as zero because xa0(0)=0 anyway
                                                       colorfactorS = 1d0          !for sQuark-Loop color already included in Loop-Function
                                                       colorfactorGl = CF 

                                                       colorfactorTU = CF                           
                                                       colorfactorQ = CF
                                                       colorfactorSH = CF           
                                                       colorfactorSV = CF

                                                       colorfactorSSVVV = CF



                        else
     
                                                       colorfactorSSVV = 1d0     
                                                       colorfactorSSVGl = 1d0 
                                                       colorfactorSSSV = 1d0 
                                                       colorfactorSSSSV = 1d0
                                                       colorfactorSSSGl = 1d0 
                                                       colorfactorQQVV = 1d0
                                                       colorfactorQQQQ = 0d0       !sQuark only contributes at O(gs^2) yet

                                                       colorfactorVS = 1d0      
                                                       colorfactorV = 1d0 
                                                       colorfactorS = 0d0          !sQuark-Loop only contributes if Looptype0=1 (only GS^2-corrections included in coulings so far)                                                      
                                                       colorfactorGl = 0d0         !no Gluino-contribution if Looptype0 not equal 1

                                                       colorfactorTU = 1d0                           
                                                       colorfactorQ = 1d0
                                                       colorfactorSH = 1d0           
                                                       colorfactorSV = 1d0

                                                       colorfactorSSVVV = 1d0

                        endif

!colorfactors appearing at tree level (including strong final states)

                        if((x_arg.eq.1).and.(y_arg.eq.1)) then                     !Gluon-Gluon final state          
                                           
                                                       colorfactorVV_SS = 4d0                                                                
                                                       colorfactorVV_TT = NC*CF**2/3d0	
                                                       colorfactorVV_TU = - 2d0/9d0                            	 
                                                       colorfactorVV_QQ = 2d0*(colorfactorVV_TT + colorfactorVV_TU) 
                                                       colorfactorVV_QTU = colorfactorVV_TT + colorfactorVV_TU
                                                       colorfactorVV_QS = 0d0                                  
                                                       colorfactorVV_STU = 2d0*cI

                        else if((x_arg.eq.1).or.(y_arg.eq.1)) then                 !Gluon-X final state  

                                                       colorfactorVH_HV = CF

                                                       colorfactorVV_SS = 0d0
                                                       colorfactorVV_TT = CF	
                                                       colorfactorVV_TU = CF                                   	 
                                                       colorfactorVV_QQ = CF                                  
                                                       colorfactorVV_QTU = CF
                                                       colorfactorVV_QS = 0d0                                  
                                                       colorfactorVV_STU = 0d0 

                        else                                                       !elw. bosons in the final state (standard case)                         

                                                       colorfactorVH_HV = 1d0  

                                                       colorfactorVV_SS = 1d0
                                                       colorfactorVV_TT = 1d0	
                                                       colorfactorVV_TU = 1d0                                   	 
                                                       colorfactorVV_QQ = 1d0                                  
                                                       colorfactorVV_QTU = 1d0
                                                       colorfactorVV_QS = 1d0                                  
                                                       colorfactorVV_STU = 1d0

                        endif

!==================================================================================================================================

        end   


