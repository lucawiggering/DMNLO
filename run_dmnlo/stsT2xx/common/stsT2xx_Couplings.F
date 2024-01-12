	subroutine stsT2xx_SetCouplings()
	implicit none

#include "stsT2xx_Model.h"
#include "stsT2xx_Kinematics.h"
#include "stsT2xx_GenCouplings.h"


      integer i, j, k, m, n, a, b, sfmeigenin1, sfmeigenout1, tt1, tt2, gen1,gen2
      double precision Krondelta
 

!*******************************
! 2-Ghost-Gluon couplings *
!*******************************
    
      GhGhG = GS
	
! ***********************************************
!2 -Squark - Gluon coupling *             
!************************************************
	
      SfSfg=0d0	
      Do i=1,2
      Do j=1,2
      Do tt1=3,4
      Do tt2=3,4
      Do gen1=1,3
      Do gen2=1,3

	SfSfg(i,tt1,gen1,j,tt2,gen2) = -GS*Krondelta(i,j)*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
 

!*********************************
! 3 gluon coupling *
!*********************************
                                      

      GGG = GS     
      
 ! ********************************
!  Fermion gluon couplings
! ********************************

      BL = 0d0
      BR = 0d0

      do i = 3, 4 
      do j = 1, 3
      do m = 3, 4
      do n = 1, 3
      BL(i,m,j,n) = -GS*Krondelta(i,m)*Krondelta(j,n)
      BR(i,m,j,n) = -GS*Krondelta(i,m)*Krondelta(j,n)
      enddo
      enddo
      enddo
      enddo

!*********************************
! 4 gluon coupling *
!*********************************
       GGGG = GS**2

!********************************************
! 2-squark-2-Gluon-couplings *
!********************************************


	SfSfgg(i,tt1,gen1,j,tt2,gen2) = 0d0

      do i=1,2
      do j=1,2
      do tt1=3,4
      do tt2=3,4
      do gen1=1,3
      do gen2=1,3
	      SfSfgg(i,tt1,gen1,j,tt2,gen2) = GS*GS*Krondelta(i,j)*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)
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
      
      

!********************************************
! 2-gluino-gluon coupling*
!********************************************
  
       GlGlG= GS

!********************************************
! 4 squarks *
!********************************************

      do sfmeigenin1=1,2
      do sfmeigenout1=1,2
      do tt1 = 1, 4
      do tt2 = 1, 4
      do gen1 = 1, 3
      do gen2 = 1, 3

      Aij(sfmeigenout1,sfmeigenin1,tt1,tt2,gen1,gen2) = (USf(sfmeigenout1,1,tt1,gen1)*USf(sfmeigenin1,1,tt1,gen1) -USf(sfmeigenout1,2,tt1,gen1)*USf(sfmeigenin1,2,tt1,gen1))*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      
	end
