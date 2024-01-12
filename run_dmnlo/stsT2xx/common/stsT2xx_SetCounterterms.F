      subroutine stsT2xx_SetCounterTerms()
      implicit none
#include "stsT2xx_Kinematics.h"
#include "stsT2xx_Model.h"
#include "stsT2xx_GenCouplings.h"
#include "stsT2xx_GenCounterterms.h"
#include "../../util/DMNLO_RenScheme.h"
      

      integer i,j,m,n,tt1,tt2,gen1,gen2
      double precision Krondelta

! Renormalization constants for testing Nf dependence
!      dZg = gs**2/(48*Pi**2)*IRdiv*4d0*4d0/2d0
!      dZSf = 0d0
!      dZC = 0d0
      

! ***********************************************
!2 -Squark -2 - Gluon counter term *             
!************************************************

      dSfSfgg=0d0
      Do i=1,2
      Do j=1,2
      Do tt1=3,4
      Do tt2=3,4
      Do gen1=1,3
      Do gen2=1,3
      dSfSfgg(i,tt1,gen1,j,tt2,gen2) = GS**2*((dZSf(j,i,tt1,gen1)+dconjg(dZSf(i,j,tt1,gen1)))/2d0+(2d0*DeltaGs+dZG)*Krondelta(i,j))*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


! ***********************************************
!2 -Squark - Gluon counter term *             
!************************************************
! i,tt1,gen1 indices of antisquark and j,tt2,gen2 indices of squark
      
	dSfSfg=0d0
      Do i=1,2
      Do j=1,2
      Do tt1=3,4
      Do tt2=3,4
      Do gen1=1,3
      Do gen2=1,3

      dSfSfg(i,tt1,gen1,j,tt2,gen2) = -GS*((dZSf(i,j,tt2,gen2)+dconjg(dZSf(j,i,tt1,gen1))+
     &dZG*Krondelta(i,j))/2d0+DeltaGS*Krondelta(i,j))*Krondelta(tt1,tt2)*Krondelta(gen1,gen2)

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
	   
!********************************************
!  3-gluon counterterm
!********************************************

      dggg=GS*(DeltaGS+3d0/2d0*dZG)
 
!********************************************
!  2-ghost-gluon counterterm
!********************************************
	   
      dGhGhG=(DeltaGs+dZC+dZG/2d0)*GS

! ********************************
!  Quark-Gluon counterterm
! ********************************

      dBL = 0d0
      dBR = 0d0

      do i = 1, 4 
      do j = 1, 4
      do m = 1, 3
      do n = 1, 3

      if((i.gt.2).and.(j.gt.2)) then

      dBL(i,j,m,n) = -GS*(DeltaGs+1/2d0*(dZG+dconjg(dZfL(j,n))+dZfL(i,m)))*Krondelta(i,j)*Krondelta(m,n)
      dBR(i,j,m,n) = -GS*(DeltaGs+1/2d0*(dZG+dconjg(dZfR(j,n))+dZfR(i,m)))*Krondelta(i,j)*Krondelta(m,n)

      endif

      enddo
      enddo
      enddo
      enddo

	end
