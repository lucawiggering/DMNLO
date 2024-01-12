!used HithcHiker-Version: 1.3

!******************************************************
!* Vector-Scalar Loop-correction to scalar propagator *      !convention(Karol/HitchHiker below. eq. 5.9)
!******************************************************

       complex*16 function QQ2xx_VS(myp22,Lptype0,Lpmeigen,Lptype,Lpgen,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS)  
                                                                                                  
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer Lptype0,Lpmeigen,Lptype,Lpgen,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1
	double precision M0Lp,M1Lp,myp22
	
* functions & subroutines

	external myBget,myAget


        M0Lp = MV(Lptype0)

        M1Lp = MSf(Lpmeigen,Lptype,Lpgen)    

        call myBget(myp22,M0Lp**2,M1Lp**2)                                            !important to keep the right order of loopfunction calls (D,C,B,A), because C calls b etc. and overwrites xb1 etc.
        call myAget(M1Lp**2)                                                                       

        g0 = VSfSf(Lptype0,Lpmeigen,LpmeigeninS,Lptype,LptypeinS,Lpgen,LpgeninS)
	g1 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen,LptypeoutS,Lptype,LpgenoutS,Lpgen)


	QQ2xx_VS =-(g0*g1*(xa0 + M0Lp**2*xb0 + 4*myp22*(xb0 + xb1)))
  
!******************************************************************************

!        print*,'QQ2xx_VS= ',QQ2xx_VS

        end
!******************************************************************************





!***************************************************
!* Vector-Scalar Loop-correction of p1 type (QQVV) *  !convention(Karol/HitchHiker below. eq. 5.109)->my conventions=correction to p1-leg
!***************************************************

       subroutine QQ2xx_QQVVp1(myp32,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,Lptype0,Lpmeigen,Lptype,Lpgen,AP,AM)  

       implicit none  
                                                                                        
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,Lptype0,Lpmeigen,Lptype,Lpgen
	complex*16 g0,g1,AP,AM
	double precision M0Lp,M1Lp,myp32
	
* functions & subroutines

	external myBget


        M0Lp = MV(Lptype0)     

        M1Lp = MSf(Lpmeigen,Lptype,Lpgen)    


        call myBget(myp32,M0Lp**2,M1Lp**2)                                                         
                                                                      
        g0 = VVSfSf(LptypeoutV,Lptype0,LpmeigenoutS,Lpmeigen,LptypeoutS,Lptype,LpgenoutS)
	g1 = VSfSf(Lptype0,Lpmeigen,LpmeigeninS,Lptype,LptypeinS,Lpgen,LpgeninS)

            
	AP = - (g0*g1*(2*xb0 + xb1))/2.
        AM =   (g0*g1*(2*xb0 + xb1))/2.

!******************************************************************************

!        print*,'AP= ',AP,'AM= ',AM

       end

!******************************************************************************




!***************************************************
!* Vector-Scalar Loop-correction of p2 type (QQVV) *   !convention(Karol/HitchHiker below. eq. 5.109)->my conventions=correction to p2-leg  
!***************************************************

       subroutine QQ2xx_QQVVp2(myp32,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,Lptype0,Lpmeigen,Lptype,Lpgen,AP,AM)  

       implicit none  
                                                                                        
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,Lptype0,Lpmeigen,Lptype,Lpgen
	complex*16 g0,g1,AP,AM
	double precision M0Lp,M1Lp,myp32
	
* functions & subroutines

	external myBget


        M0Lp = MV(Lptype0)

        M1Lp = MSf(Lpmeigen,Lptype,Lpgen)    


        call myBget(myp32,M0Lp**2,M1Lp**2)                                                         
                                                                      
        g0 = VVSfSf(Lptype0,LptypeoutV,Lpmeigen,LpmeigeninS,Lptype,LptypeinS,Lpgen)
	g1 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen,LptypeoutS,Lptype,LpgenoutS,Lpgen)

            
	AP = - (g0*g1*(2*xb0 + xb1))/2.
        AM =   (g0*g1*(2*xb0 + xb1))/2.

!******************************************************************************

!        print*,'AP= ',AP,'AM= ',AM

       end

!******************************************************************************




!***********************************************
!* Scalar-Scalar-Vector/Loop-correction (QQQQ) *      !convention(Karol-HitchHiker below eq. 5.105)
!***********************************************

       subroutine QQ2xx_VQQQQ(myp22,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,AP,AM)  

       implicit none  
                                                                                         
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1
	complex*16 g0,g1,AP,AM
	double precision M0Lp,M1Lp,colora,colorb,myp22
	
* functions & subroutines

	external myBget


        AP = (0d0,0d0)
        AM = (0d0,0d0)

        Lpgen0 = LpgeninS                                            !no generation-mixing so far
        Lpgen1 = LpgeninS

    
       do Lptype0=3,4                             
         do Lptype1=3,4


        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                 !colorconnection for up-down-mixed 4-sQuark-Vertex with Loop particles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

       do Lpmeigen0 = 1,2
         do Lpmeigen1 = 1,2

        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)    


        call myBget(myp22,M0Lp**2,M1Lp**2)                                                         

        g0 = VSfSf(LptypeoutV,Lpmeigen0,Lpmeigen1,Lptype0,Lptype1,Lpgen0,Lpgen1)


        if((LptypeinS.eq.LptypeoutS).and.(Lptype0.eq.Lptype1).and.(LptypeinS.eq.Lptype0)) then
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif


	AP = 0d0
        AM = AM + g0*g1*(xb0+2*xb1)

         enddo
       enddo

         enddo
       enddo

!******************************************************************************

!        print*,'AP= ',AP,'AM= ',AM

       end

!******************************************************************************



!***********************************************
!* Scalar-Scalar-Scalar/Loop-correction (QQQQ) *         !convention (Karol-HitchHiker below eq. 5.78)
!***********************************************

       complex*16 function QQ2xx_SQQQQ(myp22,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutH)  

       implicit none  
                                                                                      
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutH,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1
	complex*16 g0,g1,APM
	double precision M0Lp,M1Lp,colora,colorb,myp22
	
* functions & subroutines

       external myBget


       QQ2xx_SQQQQ = (0d0,0d0)

       Lpgen0 = LpgeninS                                               !no generation-mixing so far
       Lpgen1 = LpgeninS 
    
       do Lptype0=3,4                             
         do Lptype1=3,4


        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                   !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

       do Lpmeigen0 = 1,2
         do Lpmeigen1 = 1,2

        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)    


        call myBget(myp22,M0Lp**2,M1Lp**2)                                                                                                       

        g0 = HSfSf(LptypeoutH,Lpmeigen1,Lpmeigen0,Lptype1,Lptype0,Lpgen1,Lpgen0)


        if((LptypeinS.eq.LptypeoutS).and.(Lptype0.eq.Lptype1).and.(LptypeinS.eq.Lptype0)) then
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif


	APM = g0*g1*xb0

        QQ2xx_SQQQQ = QQ2xx_SQQQQ + APM

         enddo
       enddo

         enddo
       enddo

!******************************************************************************


!        print*,'QQ2xx_SQQQQ= ',QQ2xx_SQQQQ 

!******************************************************************************


       end

!******************************************************************************






!******************************************************
!* Scalar-Scalar-Scalar-Scalar/Loop-correction (QQQQ) *     !convention (Karol-HitchHiker below eq. 5.78)
!******************************************************

       complex*16 function QQ2xx_SSQQQQ(myp22,LptypeoutH1,LptypeoutH2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS)  

       implicit none  
                                                                                         
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutH1,LptypeoutH2,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1
	complex*16 g0,g1,APM
	double precision M0Lp,M1Lp,colora,colorb,myp22
	
* functions & subroutines

       external myBget


       QQ2xx_SSQQQQ = (0d0,0d0)

       Lpgen0 = LpgeninS                                             !no generation-mixing so far
       Lpgen1 = LpgeninS 
    
       do Lptype0=3,4                             
         do Lptype1=3,4

        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                 !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

       do Lpmeigen0 = 1,2
         do Lpmeigen1 = 1,2

        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)    


        call myBget(myp22,M0Lp**2,M1Lp**2)                            

        g0 = HHSfSf(LptypeoutH1,LptypeoutH2,Lpmeigen1,Lpmeigen0,Lptype1,Lptype0,Lpgen1)            

        if((LptypeinS.eq.LptypeoutS).and.(Lptype0.eq.Lptype1).and.(LptypeinS.eq.Lptype0)) then
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif

	APM = g0*g1*xb0

        QQ2xx_SSQQQQ = QQ2xx_SSQQQQ + APM

         enddo
       enddo

         enddo
       enddo

!******************************************************************************


!        print*,'QQ2xx_SSQQQQ= ',QQ2xx_SSQQQQ 

!******************************************************************************


       end

!******************************************************
!* Scalar-Scalar-Vector-Vector/Loop-correction (QQQQ) *       !convention (Karol-HitchHiker below eq. 5.105)
!******************************************************

       complex*16 function QQ2xx_VVQQQQ(myp22,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS)  

       implicit none  
                                                                                         
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV1,LptypeoutV2,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1
	complex*16 g0,g1,APM
	double precision M0Lp,M1Lp,colora,colorb,myp22
	
* functions & subroutines

       external myBget


       QQ2xx_VVQQQQ = (0d0,0d0)

       Lpgen0 = LpgeninS                                               !no generation-mixing so far
       Lpgen1 = LpgeninS 
    
       do Lptype0=3,4                             
         do Lptype1=3,4

        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                   !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

       do Lpmeigen0 = 1,2
         do Lpmeigen1 = 1,2

        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)    


        call myBget(myp22,M0Lp**2,M1Lp**2)                            

        g0 = VVSfSf(LptypeoutV1,LptypeoutV2,Lpmeigen1,Lpmeigen0,Lptype1,Lptype0,Lpgen1)            

        if((LptypeinS.eq.LptypeoutS).and.(Lptype0.eq.Lptype1).and.(LptypeinS.eq.Lptype0)) then
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,Lpmeigen1,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype1).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen1,Lptype0,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen1,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0).and.(Lptype0.eq.Lptype1)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif

	APM = g0*g1*xb0

        QQ2xx_VVQQQQ = QQ2xx_VVQQQQ + APM

         enddo
       enddo

         enddo
       enddo

!******************************************************************************


!        print*,'QQ2xx_VVQQQQ= ',QQ2xx_VVQQQQ 

!******************************************************************************


       end

!******************************************************************************





!********************************************
!* Fermion-Gluino Loop to scalar propagator *                 !convention(Karol/HitchHiker below eq. 5.11)
!********************************************

        complex*16 function QQ2xx_FGl(myp22,Lptype0,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS) 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer Lptype0,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0R,g0L,g1L,g1R,APM
	double precision myp22,M0Lp,M1Lp


* functions & subroutines

	external myBget,myAget


        if(Lptype0.eq.1) then        

        M0Lp = MGl        
 
        M1Lp = Mf(LptypeinS,LpgeninS)  

        call myBget(myp22,M0Lp**2,M1Lp**2) 
        call myAget(M1Lp**2)                                                                       
                                                                       
        g0L = GluR(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS)       !different fermion-flow-direction=>interchange PL<->PR
        g0R = GluL(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS)
        g1L = GluL(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS)
        g1R = GluR(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS)
	

	APM = - 2*((g0R*g1L + g0L*g1R)*xa0 + M0Lp*(g0R*g1L*M0Lp + g0L*g1R*M0Lp + g0L*g1L*M1Lp + g0R*g1R*M1Lp)*xb0 + (g0R*g1L + g0L*g1R)*myp22*xb1)

        else
 
        APM = (0d0,0d0)                                                     !only Gluino-correction implemented so far

        endif

!******************************************************************************

        QQ2xx_FGl = APM 

!******************************************************************************

!         print*,'QQ2xx_FGl= ',QQ2xx_FGl

        end

!******************************************************************************





!***********************************
!* Vector-Loop to scalar propagator*           !Conventions used here -> see HitchHiker below eq. 5.6
!***********************************

        complex*16 function QQ2xx_V(Lptype0,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS) 


#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer Lptype0,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1
	double precision M0Lp


* functions & subroutines

	external myAget


        M0Lp = MV(Lptype0)


        call myAget(M0Lp**2)                                                             

        g0 = VVSfSf(Lptype0,Lptype0,LpmeigenoutS,LpmeigeninS,LptypeoutS,LptypeinS,LpgenoutS)  !no generation-mixing included


	QQ2xx_V = - 2*g0*(M0Lp**2 - 2*xa0)                                                    !M0Lp=0=>xa0(0)=0 (Hitchhiker eq. 2.33, no-scale(=no mass))


!******************************************************************************

!         print*,'QQ2xx_V= ', QQ2xx_V

        end

!******************************************************************************





!***********************************
!* Scalar-Loop to scalar propagator*     !Conventions used here -> see HitchHiker above eq. 5.6 
!***********************************

        complex*16 function QQ2xx_Sfer(LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS)   !only sFermion-exchange-part proportional to GS included yet

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen0,Lptype0,Lpgen0
	complex*16 g0,g1,APM
	double precision M0Lp,colora,colorb
     
        external myAget

        QQ2xx_Sfer = (0d0,0d0)

        Lpgen0 = LpgenoutS                                           !no generation-mixing
           
      do Lpmeigen0 = 1, 2
        do Lptype0 = 3, 4                                           


        if(Lptype0.eq.LptypeinS) then

        colora = 3d0    
        colorb = 1d0                                                !colorconnection (up->down) for up-loop and down-loop

        else

        colora = 1d0    
        colorb = 3d0   

        endif


        M0Lp =  MSf(Lpmeigen0,Lptype0,Lpgen0)                     


        call myAget(M0Lp**2)                                              

!Julias Alternative

!        if(Lptype0.eq.LptypeinS) then

!        g0=4d0/3d0*sF4(Lpmeigen1,Lpmeigen0,LpmeigenoutS,Lpmeigen0,Lptype1,Lptype0,LpgenoutS)

!        else

!        g0=0d0

!        endif

        if((LptypeinS.eq.LptypeoutS).and.(LptypeinS.eq.Lptype0)) then
     
	g0 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype0,LpgenoutS)

        elseif((LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0)) then          !in/out=downtype, loop=uptype
     
	g0 =colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)

        elseif((LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0)) then          !in/out=uptype, loop=downtype
     
	g0 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype0,LpgenoutS)
   
        else
      
        g0 = 0d0
        
        endif
      
	APM = - g0*xa0

!******************************************************************************

        QQ2xx_Sfer = QQ2xx_Sfer + APM 
    
!******************************************************************************

        enddo
      enddo

!         print*,'QQ2xx_Sfer= ',QQ2xx_Sfer

	end

!******************************************************************************





!***********************************************************
!*Scalar-Scalar-Vector-Externals with Vector boson in Loop *      !follow here HitchHiker-convention below eq. 5.83
!***********************************************************

        subroutine QQ2xx_SSVV(myp12,myp22,myp32,myp13,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,AP,AM) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

 
        integer LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,AP,AM
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,myp13
	
* functions & subroutines

	external myCget, myBget

        AP = (0d0,0d0)
        AM = (0d0,0d0)

!*************************Loop-Masses******************************************
        
        M0Lp = MV(Lptype0)
        M2Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M1Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)

!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)                                                                                                       
	call myBget(myp22,M1Lp**2,M2Lp**2)
	       
	g0 = VSfSf(LptypeoutV,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1)
	g1 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen2,LptypeoutS,Lptype2,LpgenoutS,Lpgen2)
        g2 = VSfSf(Lptype0,Lpmeigen1,LpmeigeninS,Lptype1,LptypeinS,Lpgen1,LpgeninS)

	
	AP = AP + g0*g1*g2*(-4*xc00 - 2*myp12*(xc1 + xc11 + xc12) - M0Lp**2*(xc0 + xc1 + xc2) - 2*myp32*(xc12 + xc2 + xc22) - 2*myp13*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22)) 

!        print*,'AP= ',AP

        AM = AM + g0*g1*g2*(xb0 + 2*xb1 + 2*myp12*(-xc11 + xc12) + M0Lp**2*(-xc1 + xc2) + 2*myp32*(-xc12 + xc22) + 2*myp13*(-2*xc1 - xc11 + 2*xc2 + xc22))
  
!        print*,'AM= ',AM

!******************************************************************************

         end 

!******************************************************************************





!*****************************************************
!*Scalar-Scalar-Vector-Externals with Gluino in Loop *     !convention(Karol/HitchHiker below eq.5.100)->my conventions
!*****************************************************

        subroutine QQ2xx_SSVGl(myp12,myp22,myp32,myp13,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,AP,AM) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0L,g0R,g1L,g1R,g2L,g2R,AP,AM
        double precision M0Lp, M1Lp, M2Lp,myp12,myp22,myp32,myp13

        external myCget, myBget

!*************************Loop-Masses******************************************
        
        M0Lp = MGl
        M2Lp = Mf(LptypeinS,LpgeninS)
        M1Lp = Mf(LptypeoutS,LpgenoutS)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)                      
	call myBget(myp22,M1Lp**2,M2Lp**2)
	      
	g0L = BL(LptypeoutV,LptypeinS,LptypeoutS,LpgeninS,LpgenoutS)
        g0R = BR(LptypeoutV,LptypeinS,LptypeoutS,LpgeninS,LpgenoutS)
	g1L = GluL(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS)          
        g1R = GluR(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS)
        g2L = GluR(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS)               !change fermion-flow=>interchange PL<->PR
        g2R = GluL(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS) 

	
	AP =  -(g0R*(g1L*M1Lp*(g2R*M2Lp*(xc1 + xc2) + g2L*M0Lp*(xc0 + xc1 + xc2)) + g1R*(g2R*M0Lp*M2Lp*(xc0 + xc1 + xc2) + g2L*(xb0 + myp12*xc1 + myp32*xc2 + M0Lp**2*(2*xc0 + xc1 + xc2))))) - 
     &          g0L*(g1R*M1Lp*(g2L*M2Lp*(xc1 + xc2) + g2R*M0Lp*(xc0 + xc1 + xc2)) + g1L*(g2L*M0Lp*M2Lp*(xc0 + xc1 + xc2) + g2R*(xb0 + myp12*xc1 + myp32*xc2 + M0Lp**2*(2*xc0 + xc1 + xc2))))

!        print*,'AP= ',AP

        AM =    g0R*(g1L*M1Lp*(g2R*M2Lp*(-xc1 + xc2) + g2L*M0Lp*(xc0 - xc1 + xc2)) + g1R*(-(g2R*M0Lp*M2Lp*(xc0 + xc1 - xc2)) + g2L*(xb0 + 2*xb1 + myp12*xc1 - myp32*xc2 + M0Lp**2*(-xc1 + xc2)))) + 
     &          g0L*(g1R*M1Lp*(g2L*M2Lp*(-xc1 + xc2) + g2R*M0Lp*(xc0 - xc1 + xc2)) + g1L*(-(g2L*M0Lp*M2Lp*(xc0 + xc1 - xc2)) + g2R*(xb0 + 2*xb1 + myp12*xc1 - myp32*xc2 + M0Lp**2*(-xc1 + xc2))))

!        print*,'AM= ',AM

!******************************************************************************

        end

!******************************************************************************





!***********************************************************
!*Scalar-Scalar-Scalar-Externals with Vector boson in Loop *        !HitchHiker-convention above eq. 5.76 with convention(Karol)->my conventions (= M0->M1,M1->M0)
!***********************************************************

        complex*16 function QQ2xx_SSSV(myp12,myp22,myp32,myp13,LptypeoutH,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2)

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

 
        integer LptypeoutH,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lptype1,Lpgen1,Lptype2,Lpgen2
	complex*16 g0,g1,g2,APM
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,myp13

        external myCget, myBget

        QQ2xx_SSSV = (0d0,0d0)

!*************************Loop-Masses********************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)

!*******************************************************************************

	call myCget(myp12,myp22,myp32,M1Lp**2,M0Lp**2,M2Lp**2)                                              
	call myBget(myp22,M0Lp**2,M2Lp**2)
	        
        g0 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen2,LptypeoutS,Lptype2,LpgenoutS,Lpgen2)
        g1 = VSfSf(Lptype0,Lpmeigen1,LpmeigeninS,Lptype1,LptypeinS,Lpgen1,LpgeninS)
	g2 = HSfSf(LptypeoutH,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        
	
	APM = g0*g1*g2*(xb0 + M1Lp**2*xc0 + myp12*(xc0 - 2*xc1) + 2*myp32*xc2 - 2*myp13*(xc0 - xc1 + xc2))  

!******************************************************************************
       
        QQ2xx_SSSV = QQ2xx_SSSV + APM

!******************************************************************************

!         print*,'QQ2xx_SSSV= ',QQ2xx_SSSV

         end   

!******************************************************************************                                                                





!******************************************************************
!*Scalar-Scalar-Scalar-Scalar-Externals with Vector boson in Loop *       !HitchHiker-convention above eq. 5.76 with convention(Karol)->my conventions (= M0->M1,M1->M0)
!******************************************************************

        complex*16 function QQ2xx_SSSSV(myp12,myp22,myp32,myp13,LptypeoutH1,LptypeoutH2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2)

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

 
        integer LptypeoutH1,LptypeoutH2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype0,Lptype1,Lpgen1,Lptype2,Lpgen2
	complex*16 g0,g1,g2,APM
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,myp13

        external myCget, myBget

        QQ2xx_SSSSV = (0d0,0d0)

!*************************Loop-Masses********************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)

!*******************************************************************************

	call myCget(myp12,myp22,myp32,M1Lp**2,M0Lp**2,M2Lp**2)                                              
	call myBget(myp22,M0Lp**2,M2Lp**2)
	        
        g0 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen2,LptypeoutS,Lptype2,LpgenoutS,Lpgen2)
        g1 = VSfSf(Lptype0,Lpmeigen1,LpmeigeninS,Lptype1,LptypeinS,Lpgen1,LpgeninS)
	g2 = HHSfSf(LptypeoutH1,LptypeoutH2,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1)
        
	
	APM = g0*g1*g2*(xb0 + M1Lp**2*xc0 + myp12*(xc0 - 2*xc1) + 2*myp32*xc2 - 2*myp13*(xc0 - xc1 + xc2))  

!******************************************************************************
       
        QQ2xx_SSSSV = QQ2xx_SSSSV + APM

!******************************************************************************

!         print*,'QQ2xx_SSSSV= ',QQ2xx_SSSSV

         end   

!******************************************************************************   





!******************************************************************
!*Scalar-Scalar-Vector-Vector-Externals with Vector boson in Loop *    !Conventions as in HitchHiker below eq. 5.83
!******************************************************************

       complex*16 function QQ2xx_VVQQVV(myp12,myp22,myp32,myp13,Lptype0,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2)  

       implicit none  
                                                                                         
                                                                                                                   
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer Lptype0,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutV1,LptypeoutV2,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,APM
	double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,myp13
	
* functions & subroutines

       external myBget


       QQ2xx_VVQQVV = (0d0,0d0)

       M0Lp = MV(Lptype0)

    
       M1Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)                                    
       M2Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)      

       call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)
       call myBget(myp22,M1Lp**2,M2Lp**2)                            

       g0 = VVSfSf(LptypeoutV1,LptypeoutV2,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2)
       g1 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen2,LptypeoutS,Lptype2,LpgenoutS,Lpgen2)
       g2 = VSfSf(Lptype0,Lpmeigen1,LpmeigeninS,Lptype1,LptypeinS,Lpgen1,LpgeninS)            

       APM = g0*g1*g2*(xb0 + M0Lp**2*xc0 + 2*(myp12*xc1 + myp32*xc2 + myp13*(2*xc0 + xc1 + xc2)))

       QQ2xx_VVQQVV = QQ2xx_VVQQVV + APM

!******************************************************************************


!        print*,'QQ2xx_VVQQVV= ',QQ2xx_VVQQVV 

!******************************************************************************


       end

!******************************************************************************





!*************************************************************************
!*Scalar-Scalar-Scalar-Externals with Gluino in Loop (Gluino,Neutralino) *            !convention(Karol/HitchHiker below eq. 5.71)
!*************************************************************************

        complex*16 function QQ2xx_SSSGl(myp12,myp22,myp32,myp13,LptypeoutH,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS) 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,LptypeoutH
	complex*16 g0L,g0R,g1L,g1R,g2L,g2R,APM
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,myp13

        external myCget, myBget
 
!*************************Loop-Masses******************************************
        
        M0Lp = MGl
        M2Lp = Mf(LptypeinS,LpgeninS)
        M1Lp = Mf(LptypeoutS,LpgenoutS)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)                  
	call myBget(myp22,M1Lp**2,M2Lp**2)
        
	g0L = DL(LptypeoutH,LptypeinS,LptypeoutS,LpgeninS,LpgenoutS)
        g0R = DR(LptypeoutH,LptypeinS,LptypeoutS,LpgeninS,LpgenoutS)
	g1L = GluL(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS)           
        g1R = GluR(LpmeigenoutS,LptypeoutS,LptypeoutS,LpgenoutS,LpgenoutS) 
        g2L = GluR(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS)           !change fermion-flow=>interchange PL<->PR
        g2R = GluL(LpmeigeninS,LptypeinS,LptypeinS,LpgeninS,LpgeninS)

	
	APM =  2*(M0Lp*((g0L*g1L*g2L + g0R*g1R*g2R)*M1Lp*M2Lp + (g0R*g1L*g2L + g0L*g1R*g2R)*myp13)*xc0 + 
     &    (g0L*(g1R*g2R*M0Lp + g1L*g2R*M1Lp + g1R*g2L*M2Lp) + g0R*(g1L*g2L*M0Lp + g1R*g2L*M1Lp + g1L*g2R*M2Lp))*(xb0 + M0Lp**2*xc0) + 
     &    (g0L*g1R*(g2R*M0Lp + g2L*M2Lp) + g0R*g1L*(g2L*M0Lp + g2R*M2Lp))*(myp12*xc1 + myp13*xc2) + (g0L*g2R*(g1R*M0Lp + g1L*M1Lp) + g0R*g2L*(g1L*M0Lp + g1R*M1Lp))*(myp13*xc1 + myp32*xc2))

!******************************************************************************

        QQ2xx_SSSGl = APM

!******************************************************************************

!       print *, "QQ2xx_SSSGl=", QQ2xx_SSSGl 

	end

!*****************************************************************************





!***************************
!*Box-Triangle (Typ1,HV(t))*        HitchHiker convention below eq.5.75 with p1-attached coupling replaced by Scalar-Scalar-Vector-Vector-coupling
!***************************

        subroutine QQ2xx_BT1HV(tu_switch,myp12,myp22,myp32,Lptype0,LptypeoutH,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,A1,A3) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,Lptype0,LptypeoutH,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,A1,A3
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32

        external myCget
        
        A1 = (0d0,0d0)
        A3 = (0d0,0d0)

        do Lpmeigen1 = 1,2

!*************************Loop-Masses******************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)  

        if (tu_switch.eq.0) then                                !chose t- or u-channel

	g0 = HSfSft(LptypeoutH,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1)
        g1 = VVSfSfu(LptypeoutV,Lptype0,Lpmeigen1,LpmeigenoutS,Lptype1,LptypeoutS,LpgenoutS)

        else

	g0 = HSfSfu(LptypeoutH,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1)
        g1 = VVSfSft(LptypeoutV,Lptype0,Lpmeigen1,LpmeigenoutS,Lptype1,LptypeoutS,LpgenoutS)

        endif

	g2 = VSfSf(Lptype0,Lpmeigen2,LpmeigeninS,Lptype2,LptypeinS,Lpgen2,LpgeninS)
       
 	
	A1 = A1 - g0*g1*g2*xc1                        !p1-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A3 = A3 - g0*g1*g2*(2*xc0 + xc2)              !p3-coefficient->see mathematica-file
!        print*,'A3= ',A3

!******************************************************************************

        enddo

        end

!******************************************************************************





!***************************
!*Box-Triangle (Typ2,VH(t))*     HitchHiker convention below eq.5.75 with coupling-name g1<->g2 interchanged and p3-attached coupling replaced by Scalar-Scalar-Vector-Vector-coupling
!***************************

        subroutine QQ2xx_BT2VH(tu_switch,myp12,myp22,myp32,Lptype0,LptypeoutH,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lptype2,Lpgen2,A1,A3) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,Lptype0,LptypeoutH,LptypeoutV,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,A1,A3
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32

        external myCget
        
        A1 = (0d0,0d0)
        A3 = (0d0,0d0)

        do Lpmeigen2 = 1,2

!*************************Loop-Masses******************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)  

        if (tu_switch.eq.0) then                                !chose t- or u-channel

	g0 = HSfSfu(LptypeoutH,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1)
        g1 = VVSfSft(LptypeoutV,Lptype0,LpmeigeninS,Lpmeigen2,LptypeinS,Lptype2,LpgeninS)

        else

	g0 = HSfSft(LptypeoutH,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1)
        g1 = VVSfSfu(LptypeoutV,Lptype0,LpmeigeninS,Lpmeigen2,LptypeinS,Lptype2,LpgeninS)

        endif

	g2 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen1,LptypeoutS,Lptype1,LpgenoutS,Lpgen1)
        	

	A1 = A1 - (g0*g1*g2*(2*xc0 + xc1))      !p1-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A3 = A3 - (g0*g1*g2*xc2)                !p3-coefficient->see mathematica-file
!        print*,'A3= ',A3

!******************************************************************************

        enddo

        end

!******************************************************************************





!************************
!*Box-Triangle (Typ1,VV)*       HitchHiker convention below eq.5.83 with p1-attached coupling replaced by Scalar-Scalar-Vector-Vector-coupling
!************************

        subroutine QQ2xx_BT1VV(tu_switch,myp12,myp22,myp32,Lptype0,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,A0,A1,A2,A3,A4) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,Lptype0,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,A0,A1,A2,A3,A4
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32

        external myCget
        
        A0 = (0d0,0d0)
        A1 = (0d0,0d0)
        A2 = (0d0,0d0)
        A3 = (0d0,0d0)
        A4 = (0d0,0d0)


        do Lpmeigen1 = 1,2

!*************************Loop-Masses******************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)  

        if (tu_switch.eq.0) then                                !chose t- or u-channel

	g0 = VSfSft(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g1 = VVSfSfu(LptypeoutV2,Lptype0,Lpmeigen1,LpmeigenoutS,Lptype1,LptypeoutS,LpgenoutS)

        else

       	g0 = VSfSfu(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g1 = VVSfSft(LptypeoutV2,Lptype0,Lpmeigen1,LpmeigenoutS,Lptype1,LptypeoutS,LpgenoutS)

        endif

	g2 = VSfSf(Lptype0,Lpmeigen2,LpmeigeninS,Lptype2,LptypeinS,Lpgen2,LpgeninS)
        	

        A0 = A0 + 2*g0*g1*g2*xc00                               !g_mu,nu-coefficient->see mathematica-file
!        print*,'A0= ',A0

        A1 = A1 + (g0*g1*g2*(xc1 + 2*xc11))                     !p1_mu*p1_nu-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A2 = A2 + (g0*g1*g2*(2*xc0 + 4*xc1 + 2*xc12 + xc2))     !p1_mu*p3_nu-coefficient->see mathematica-file   
!        print*,'A2= ',A2

        A3 = A3 + (g0*g1*g2*(xc1 + 2*xc12))                     !p3_mu*p1_nu-coefficient->see mathematica-file
!        print*,'A3= ',A3

        A4 = A4 + (g0*g1*g2*(2*xc0 + 5*xc2 + 2*xc22))           !p3_mu*p3_nu-coefficient->see mathematica-file
!        print*,'A4= ',A4



!******************************************************************************

        enddo

        end

!******************************************************************************





!************************
!*Box-Triangle (Typ2,VV)*              HitchHiker convention below eq.5.83 with coupling-name g1<->g2 interchanged and p3-attached coupling replaced by Scalar-Scalar-Vector-Vector-coupling 
!************************

        subroutine QQ2xx_BT2VV(tu_switch,myp12,myp22,myp32,Lptype0,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lptype2,Lpgen2,A0,A1,A2,A3,A4) 

        implicit none 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,Lptype0,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2
	complex*16 g0,g1,g2,A0,A1,A2,A3,A4
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32

        external myCget
                
        A0 = (0d0,0d0)
        A1 = (0d0,0d0)
        A2 = (0d0,0d0)
        A3 = (0d0,0d0)
        A4 = (0d0,0d0)

        do Lpmeigen2 = 1,2

!*************************Loop-Masses******************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)          

        if (tu_switch.eq.0) then                                         !chose t- or u-channel

	g0 = VSfSfu(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g1 = VVSfSft(LptypeoutV2,Lptype0,LpmeigeninS,Lpmeigen2,LptypeinS,Lptype2,LpgeninS)

        else

       	g0 = VSfSft(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g1 = VVSfSfu(LptypeoutV2,Lptype0,LpmeigeninS,Lpmeigen2,LptypeinS,Lptype2,LpgeninS)

        endif

	g2 = VSfSf(Lptype0,LpmeigenoutS,Lpmeigen1,LptypeoutS,Lptype1,LpgenoutS,Lpgen1)
        	

        A0 = A0 + 2*g0*g1*g2*xc00                               !g_mu,nu-coefficient->see mathematica-file
!        print*,'A0= ',A0

        A1 = A1 + g0*g1*g2*(2*xc0 + 5*xc1 + 2*xc11)             !p1_mu*p1_nu-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A2 = A2 + g0*g1*g2*(2*xc12 + xc2)                       !p1_mu*p3_nu-coefficient->see mathematica-file
!        print*,'A2= ',A2

        A3 = A3 + g0*g1*g2*(2*xc0 + xc1 + 2*xc12 + 4*xc2)       !p3_mu*p1_nu-coefficient->see mathematica-file
!        print*,'A3= ',A3

        A4 = A4 + g0*g1*g2*(xc2 + 2*xc22)                       !p3_mu*p3_nu-coefficient->see mathematica-file
!        print*,'A4= ',A4

!******************************************************************************

        enddo

        end

!******************************************************************************





!***********************
!*box-selfenergy-hybrid*  usual t/u-channel-conventions with two Scalar-Scalar-Vector-Vector-couplings
!***********************

        complex*16 function QQ2xx_BSVV(tu_switch,myp12,Lptype0,LptypeoutV1,LptypeoutV2,Lptype,Lpgen,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS) 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,Lptype0,LptypeoutV1,LptypeoutV2,Lpmeigen,Lptype,Lpgen,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1,A0
        double precision M0Lp,M1Lp,myp12

        external myBget
 

        A0=(0d0,0d0)


        do Lpmeigen = 1,2

!*************************Loop-Masses********************************************
        
        M0Lp = MV(Lptype0)
        M1Lp = MSf(Lpmeigen,Lptype,Lpgen)
 
!******************************************************************************

	call myBget(myp12,M0Lp**2,M1Lp**2)                   

        if (tu_switch.eq.0) then                             !chose t- or u-channel   

	g0 = VVSfSft(Lptype0,LptypeoutV1,LpmeigeninS,Lpmeigen,LptypeinS,Lptype,LpgeninS)
        g1 = VVSfSfu(Lptype0,LptypeoutV2,Lpmeigen,LpmeigenoutS,Lptype,LptypeoutS,Lpgen)

        else

	g0 = VVSfSfu(Lptype0,LptypeoutV1,LpmeigeninS,Lpmeigen,LptypeinS,Lptype,LpgeninS)
        g1 = VVSfSft(Lptype0,LptypeoutV2,Lpmeigen,LpmeigenoutS,Lptype,LptypeoutS,Lpgen)

        endif
      
	
	A0 =  A0 - g0*g1*xb0

!******************************************************************************

        enddo


        QQ2xx_bsVV = A0

!******************************************************************************

!       print *, "QQ2xx_bsVV=", QQ2xx_bsVV 

	end

 !*****************************************************************************






!*****************************************************************
!*box-triangle-with 4-sQuark-Vertex for Scalar-Scalar-final-state*   HitchHiker convention below eq.5.75 with p1-attached coupling replaced by 4-sQuark-coupling
!*****************************************************************

        complex*16 function QQ2xx_BTsf4SS(tu_switch,myp12,myp22,myp32,LptypeoutH1,LptypeoutH2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS) 

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,LptypeoutH1,LptypeoutH2,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1,g2,A0
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,colora,colorb

        external myCget
 

        A0=(0d0,0d0)
       
        Lpgen0=LpgeninS                                                   !no generation-mixing
        Lpgen1=LpgenoutS
        Lpgen2=LpgeninS

        do Lptype0=3,4                                               
         do Lptype1=3,4

        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                     !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

          do Lptype2=3,4   
      
            do Lpmeigen0 = 1,2
             do Lpmeigen1 = 1,2
               do Lpmeigen2 = 1,2

!*************************Loop-Masses********************************************
        
        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!********************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)           

        if (tu_switch.eq.0) then                                         !chose t- or u-channel  

	g0 = HSfSfu(LptypeoutH1,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1) 
        g2 = HSfSft(LptypeoutH2,Lpmeigen0,Lpmeigen2,Lptype0,Lptype2,Lpgen0,Lpgen2)

        else

	g0 = HSfSft(LptypeoutH1,Lpmeigen2,Lpmeigen1,Lptype2,Lptype1,Lpgen2,Lpgen1) 
        g2 = HSfSfu(LptypeoutH2,Lpmeigen0,Lpmeigen2,Lptype0,Lptype2,Lpgen0,Lpgen2)

        endif


        if((LptypeinS.eq.LptypeoutS).and.(Lptype1.eq.Lptype0).and.(LptypeinS.eq.Lptype1)) then                                                   !Lptype0 runs out of the 4-squark-vertex, Lptype1 goes in.
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif

        	
	A0 =  A0 - g0*g1*g2*xc0

!******************************************************************************

             enddo
            enddo
           enddo
          enddo
         enddo
        enddo



        QQ2xx_BTsf4SS = A0

!******************************************************************************

!       print *, "QQ2xx_BTsf4SS=", QQ2xx_BTsf4SS

	end

 !*****************************************************************************






!*****************************************************************
!*box-triangle-with 4-sQuark-Vertex for Vector-Vector-final-state*      HitchHiker convention below eq.5.83 with p3-leg taken as a Vector and p1-attached coupling replaced by 4-sQuark-coupling
!*****************************************************************

        subroutine QQ2xx_BTsf4VV(tu_switch,myp12,myp22,myp32,LptypeoutV1,LptypeoutV2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,A0,A1,A2,A3,A4) 

        implicit none  

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,LptypeoutV1,LptypeoutV2,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1,g2,A0,A1,A2,A3,A4
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,colora,colorb

        external myCget
 

        A0 = (0d0,0d0)
        A1 = (0d0,0d0)
        A2 = (0d0,0d0)
        A3 = (0d0,0d0)
        A4 = (0d0,0d0)
       
        Lpgen0=LpgeninS                                                  !no generation-mixing
        Lpgen1=LpgenoutS
        Lpgen2=LpgeninS

        do Lptype0=3,4                                               
         do Lptype1=3,4

        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                    !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

          do Lptype2=3,4   
      
            do Lpmeigen0 = 1,2
             do Lpmeigen1 = 1,2
               do Lpmeigen2 = 1,2

!*************************Loop-Masses********************************************
        
        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)           

        if (tu_switch.eq.0) then                                         !chose t- or u-channel 

	g0 = VSfSfu(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = VSfSft(LptypeoutV2,Lpmeigen2,Lpmeigen0,Lptype2,Lptype0,Lpgen2,Lpgen0)

        else

	g0 = VSfSft(LptypeoutV1,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = VSfSfu(LptypeoutV2,Lpmeigen2,Lpmeigen0,Lptype2,Lptype0,Lpgen2,Lpgen0)

        endif        


        if((LptypeinS.eq.LptypeoutS).and.(Lptype1.eq.Lptype0).and.(LptypeinS.eq.Lptype1)) then                                                   !Lptype0 runs out of the 4-squark-vertex, Lptype1 goes in.
     
        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif

	
        A0 = A0 - 4*g0*g1*g2*xc00                               !g_mu,nu-coefficient->see mathematica-file
!        print*,'A0= ',A0

        A1 = A1 - 2*g0*g1*g2*(xc1 + 2*xc11)                     !p1_mu*p1_nu-coefficient->see mathematica-file                     
!        print*,'A1= ',A1

        A2 = A2 - (g0*g1*g2*(xc0 + 2*(xc1 + 2*xc12 + xc2)))     !p1_mu*p3_nu-coefficient->see mathematica-file
!        print*,'A2= ',A2

        A3 = A3 - 2*g0*g1*g2*(xc1 + 2*xc12)                     !p3_mu*p1_nu-coefficient->see mathematica-file
!        print*,'A3= ',A3

        A4 = A4 - (g0*g1*g2*(xc0 + 4*(xc2 + xc22)))             !p3_mu*p3_nu-coefficient->see mathematica-file
!        print*,'A4= ',A4

!******************************************************************************

             enddo
            enddo
           enddo
          enddo
         enddo
        enddo


	end

 !*****************************************************************************





!***********************************************************************
!*box-triangle-with 4-sQuark-Vertex for Vector-Higgs-final-state (Typ1)*         HitchHiker convention below eq.5.83 with p1-attached coupling replaced by 4-sQuark-coupling
!***********************************************************************

        subroutine QQ2xx_BTsf4VH1(tu_switch,myp12,myp22,myp32,LptypeoutV,LptypeoutH,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,A1,A3) 

        implicit none  

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,LptypeoutV,LptypeoutH,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1,g2,A1,A3
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,colora,colorb

        external myCget
 

        A1 = (0d0,0d0)
        A3 = (0d0,0d0)
       
        Lpgen0=LpgeninS                             !no generation-mixing
        Lpgen1=LpgenoutS
        Lpgen2=LpgeninS

        do Lptype0=3,4                                               
         do Lptype1=3,4

        if((Lptype0.eq.Lptype1).and.(LptypeinS.ne.Lptype1)) then

        colora = 1d0 
        colorb = 3d0                                                 !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

          do Lptype2=3,4   
      
            do Lpmeigen0 = 1,2
             do Lpmeigen1 = 1,2
               do Lpmeigen2 = 1,2

!*************************Loop-Masses********************************************
        
        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)           

        if (tu_switch.eq.0) then                                         !chose t- or u-channel  

	g0 = VSfSfu(LptypeoutV,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = HSfSft(LptypeoutH,Lpmeigen0,Lpmeigen2,Lptype0,Lptype2,Lpgen0,Lpgen2)

        else

	g0 = VSfSft(LptypeoutV,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = HSfSfu(LptypeoutH,Lpmeigen0,Lpmeigen2,Lptype0,Lptype2,Lpgen0,Lpgen2)

        endif 


        if((LptypeinS.eq.LptypeoutS).and.(Lptype1.eq.Lptype0).and.(LptypeinS.eq.Lptype1)) then                                                   !Lptype0 comes out of the 4-squark-vertex, Lptype1 goes in.
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,Lpmeigen0,LpmeigenoutS,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype0).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype1)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,LpmeigenoutS,Lpmeigen0,Lptype1,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen1,LpmeigenoutS,Lpmeigen0,LptypeinS,Lptype1,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype1).and.(Lptype1.eq.Lptype0)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen1,LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lptype1,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif


        A1 = A1 + g0*g1*g2*(xc0 + 2*xc1)           !p1-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A3 = A3 + g0*g1*g2*(xc0 + 2*xc2)           !p3-coefficient->see mathematica-file
!        print*,'A3= ',A3


!******************************************************************************

             enddo
            enddo
           enddo
          enddo
         enddo
        enddo


	end

 !*****************************************************************************





!***********************************************************************
!*box-triangle-with 4-sQuark-Vertex for Vector-Higgs-final-state (Typ2)*   HitchHiker convention below eq.5.83 with coupling names g1<->g2 interchanged and p3-attached coupling replaced by 4-sQuark-coupling
!***********************************************************************

        subroutine QQ2xx_BTsf4VH2(tu_switch,myp12,myp22,myp32,LptypeoutV,LptypeoutH,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS,A1,A3) 

        implicit none  

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer tu_switch,LptypeoutV,LptypeoutH,Lpmeigen0,Lptype0,Lpgen0,Lpmeigen1,Lptype1,Lpgen1,Lpmeigen2,Lptype2,Lpgen2,LpmeigeninS,LptypeinS,LpgeninS,LpmeigenoutS,LptypeoutS,LpgenoutS
	complex*16 g0,g1,g2,A1,A3
        double precision M0Lp,M1Lp,M2Lp,myp12,myp22,myp32,colora,colorb

        external myCget
 

        A1 = (0d0,0d0)
        A3 = (0d0,0d0)
       
        Lpgen0=LpgeninS                                                  !no generation-mixing
        Lpgen1=LpgenoutS
        Lpgen2=LpgeninS

        do Lptype0=3,4                                               
         do Lptype2=3,4

        if((Lptype2.eq.Lptype0).and.(LptypeinS.ne.Lptype0)) then

        colora = 1d0 
        colorb = 3d0                                                    !colorconnection for up-down-mixed 4-sQuark-Vertex with Loopparticles of same type and color

        else

        colora = 3d0                                                 
        colorb = 1d0    

        endif

          do Lptype1=3,4   
      
            do Lpmeigen0 = 1,2
             do Lpmeigen1 = 1,2
               do Lpmeigen2 = 1,2

!*************************Loop-Masses********************************************
        
        M0Lp = MSf(Lpmeigen0,Lptype0,Lpgen0)
        M1Lp = MSf(Lpmeigen1,Lptype1,Lpgen1)
        M2Lp = MSf(Lpmeigen2,Lptype2,Lpgen2)
 
!******************************************************************************

	call myCget(myp12,myp22,myp32,M0Lp**2,M1Lp**2,M2Lp**2)           

        if (tu_switch.eq.0) then                                         !chose t- or u-channel  

	g0 = VSfSft(LptypeoutV,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = HSfSfu(LptypeoutH,Lpmeigen1,Lpmeigen0,Lptype1,Lptype0,Lpgen1,Lpgen0)

        else

	g0 = VSfSfu(LptypeoutV,Lpmeigen1,Lpmeigen2,Lptype1,Lptype2,Lpgen1,Lpgen2)
        g2 = HSfSft(LptypeoutH,Lpmeigen1,Lpmeigen0,Lptype1,Lptype0,Lpgen1,Lpgen0)

        endif 

       
        if((LptypeinS.eq.LptypeoutS).and.(Lptype0.eq.Lptype2).and.(LptypeinS.eq.Lptype0)) then                                                   !Lptype2 comes out of the 4-squark-vertex, Lptype0 goes in.
     
        g1 =colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen2,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen2,LptypeinS,Lptype0,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.Lptype2).and.(4.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming up, outgoing down

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,Lpmeigen2,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,Lpmeigen2,LpmeigenoutS,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.Lptype2).and.(3.eq.LptypeoutS).and.(LptypeoutS.eq.Lptype0)) then                               !incoming down, outgoing up

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen2,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,LpmeigenoutS,Lpmeigen2,Lptype0,LptypeinS,LpgenoutS)

        elseif((3.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(4.eq.Lptype0).and.(Lptype0.eq.Lptype2)) then                                  !incoming up, outgoing up

        g1 = colora*sF4a(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen2,LptypeinS,Lptype0,LpgenoutS)+colorb*sF4b(LpmeigeninS,Lpmeigen0,LpmeigenoutS,Lpmeigen2,LptypeinS,Lptype0,LpgenoutS)

        elseif((4.eq.LptypeinS).and.(LptypeinS.eq.LptypeoutS).and.(3.eq.Lptype0).and.(Lptype0.eq.Lptype2)) then                                  !incoming down, outgoing down

        g1 = colora*sF4a(Lpmeigen0,LpmeigeninS,Lpmeigen2,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)+colorb*sF4b(Lpmeigen0,LpmeigeninS,Lpmeigen2,LpmeigenoutS,Lptype0,LptypeinS,LpgenoutS)
    
        else
      
        g1 = 0d0
        
        endif


        A1 = A1 + g0*g1*g2*(xc0 + 2*xc1)                        !p1-coefficient->see mathematica-file
!        print*,'A1= ',A1

        A3 = A3 + g0*g1*g2*(xc0 + 2*xc2)                        !p3-coefficient->see mathematica-file
!        print*,'A3= ',A3


!******************************************************************************

             enddo
            enddo
           enddo
          enddo
         enddo
        enddo


	end

 !*****************************************************************************