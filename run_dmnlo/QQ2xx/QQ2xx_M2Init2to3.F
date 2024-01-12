!**********************************************************************************************************************************

        subroutine QQ2xx_Init2to3()

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_Init2to3.h"
#include "QQ2xx_GenCouplings.h"

	integer i

!==================================================================================================================================
!                                           ----Propagatoren----
!==================================================================================================================================

!only for initial-state-radiation-> s-Channel: sqrtS^2 -> (p1+p2-k3)^2 = sqrtS**2*(-1d0 + muu5**2 + x1 + x2) 

	DenSmV23(1)  =  sqrtS**2*(-1d0 + muu5**2 + x1 + x2)                                                     !Gluon

	DenSmV23(2)  =  sqrtS**2*(-1d0 + muu5**2 + x1 + x2)                                                     !Photon

	DenSmV23(3)  =  (sqrtS**2*(-1d0 + muu5**2 + x1 + x2) - MZ2) + cI*MZ*wZ                                  !Z-Boson

	DenSmV23(4)  =  (sqrtS**2*(-1d0 + muu5**2 + x1 + x2) - MW2) + cI*MW*wW                                  !W-Boson

!=============================----Index:1=h,2=H,3=P,4=G0,5=H+,6=G+----=============================================================
 
        do i = sumstartH, sumendH                                                    
                                                             
	  DenSmH23(i) = (sqrtS**2*(-1d0 + muu5**2 + x1 + x2) - MHk(i)**2) + cI*MHk(i)*wHk(i)                    !Higgs-Boson

	enddo

!==================================================================================================================================    
!                                            ----Mandelstam- and other variables---- 
!==================================================================================================================================

                                                       if(ztype.eq.1) then
                                
                                                       colorfactor23 = CF

                                                       else
 
                                                       colorfactor23 = 1d0

                                                       endif             

                                    
                               V1rad = VSfSf(ztype,sfmeigen1,sfmeigen1,sftype1,sftype1,sfgen1,sfgen1)           !p1-radiation
                               V2rad = VSfSf(ztype,sfmeigen2,sfmeigen2,sftype2,sftype2,sfgen2,sfgen2)           !p2-radiation   

                               V1radc = dconjg(VSfSf(ztype,sfmeigen1,sfmeigen1,sftype1,sftype1,sfgen1,sfgen1))  
                               V2radc = dconjg(VSfSf(ztype,sfmeigen2,sfmeigen2,sftype2,sftype2,sfgen2,sfgen2))        

                           
!==================================================================================================================================    
!                                                     ----on shell check----                                                                          
!================================================================================================================================== 

                               do i=1,4
!  
!                                    if(dabs(dreal(DenSmV23(i))).lt.1d0) then           !set on shell contributions already included in 2->2 part to zero
! 
!                                       onshell_V(i) = 0d0
! 
!                                    else
! 
                                      onshell_V(i) = 1d0
! 
!                                    endif
! 
                               enddo
! 
! 
                               do i=sumstartH, sumendH
!  
!                                   if(dabs(dreal(DenSmH23(i))).lt.1d0) then           !set on shell contributions already included in 2->2 part to zero
! 
!                                     onshell_H(i) = 0d0
! 
!                                   else
! 
                                     onshell_H(i) = 1d0
! 
!                                   endif
! 


                               enddo

!==================================================================================================================================

        end   


