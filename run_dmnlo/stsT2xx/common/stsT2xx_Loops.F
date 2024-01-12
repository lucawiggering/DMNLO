! The subroutine "stsT2xx_calc_partial_amplitudes_ASqSqG" computes the form factors Am and Ap for the squark-squark-gluon vertex
! which are defined in the MSc Thesis by Luca. All external particles are treated as in general off-shell.
      subroutine stsT2xx_calc_partial_amplitudes_ASqSqG(Am,Ap,sf1,tt1,gen1,sf2,tt2,gen2,p1p1,p2p2,p1Pp22)
! sf1,tt1,gen1 are sfermion index, type and generation of squark with momentum p1 (scalar product p1^2=p1p1).
! sf2,tt2,gen2 are sfermion index, type and generation of antisquark  with momentum p2 (scalar product p2^2=p2p2).
! The gluon carries the squared momentum (p1+p2)^2=p1Pp22.
! Only Am contains the UV divergences due to its proportionality to the tree level vertex. 
      implicit none
#include "stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"
       
      double complex Am,Ap,mixing
      integer sf1,tt1,gen1,sf2,tt2,gen2
      double precision color,p1p1,p2p2,p1Pp22,p1p2,M2,MQ,MQ2,Krondelta

      M2=MSf2(sf1,tt1,gen1)

      if (gen1.eq.3) then 
      MQ=Mf(tt1,gen1)
      MQ2=Mf2(tt1,gen1)
      else
      MQ=0d0
      MQ2=0d0
      endif

      Am=0d0
      Ap=0d0
 
      p1p2=(p1Pp22-p1p1-p2p2)/2d0
      M2=MSf2(sf1,tt1,gen1)
      mixing=USf(sf1,2,tt1,gen1)*USf(sf2,1,tt1,gen1)+USf(sf1,1,tt1,gen1)*USf(sf2,2,tt1,gen1)

      if (tt1.eq.tt2 .and. gen1.eq.gen2) then

      if (sf1.eq.sf2) then
      
! Gluon-squark loop at p1 (1)
      call myBget(p1p1,M2,0d0)

      color=(Nc**2-2d0)/(2d0*Nc)
      Ap=Ap+color*(GS**3)/(32d0*Pi**2)*(xb0-xb1) 
      Am=Am+color*(GS**3)/(32d0*Pi**2)*(xb0-xb1) 
 
! Gluon-squark loop at p2 (2)
      call myBget(p2p2,M2,0d0)

      color=(Nc**2-2d0)/(2d0*Nc)
      Ap=Ap+color*(GS**3)/(32d0*Pi**2)*(xb1-xb0) 
      Am=Am-color*(GS**3)/(32d0*Pi**2)*(xb1-xb0)

! Squark loop (3)
!     yields zero

! Gluon loop (4)
!     yields zero
 
! Gluon exchange (5) (Gluon-Squark-Squark)
      call myCget(p1p1,p1Pp22,p2p2,0d0,M2,M2)
      call myBget(p1Pp22,M2,M2)

      color=-1d0/(Nc*2d0)
      Ap=Ap+color*(GS**3*(xb0 + 2*xb1 - 2*p1p1*(xc11 - xc12) + 2*p2p2*(xc22-xc12)+ 2*p1p2*(2*xc1 + xc11 - 2*xc2 - xc22) ))/(16d0*Pi**2)
      Am=Am-color*(GS**3*(2*xc00 + p2p2*(xc12+ xc2 +xc22) + p1p1*(xc1 + xc11+ xc12) - p1p2*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22)))/(8d0*Pi**2)

! Squark exchange (6) (Squark-Gluon-Gluon)
      call myCget(p1p1,p1Pp22,p2p2,M2,0d0,0d0) 
      call myBget(p1Pp22,0d0,0d0)
      color=Nc/2d0
      Ap=Ap+color*(GS**3*(p1p1*(xc0 - 3*xc1 + 2*xc11 - 2*xc12 - xc2) + p2p2*(-xc0 + xc1 + 2*xc12 + 3*xc2 - 2*xc22) + 2*p1p2*(xc1 - xc11 - xc2 + xc22)))/(32d0*Pi**2)
      Am=Am+color*(GS**3*(-4*xb0 + 4*xc00 + p2p2*(-xc0 + xc1 + 2*xc12 + 3*xc2 + 2*xc22) + p1p1*(-xc0 + 3*xc1 + 2*xc11 + 2*xc12 + xc2) - 
     &     2*p1p2*(xc0 - 2*xc1 + xc11 + 2*xc12 - 2*xc2 + xc22) - 4*xc0*M2))/(32d0*Pi**2)

      endif

! Gluino exchange (7) (Gluino-Quark-Quark) 
      call myCget(p1p1,p1Pp22,p2p2,MGl2,MQ2,MQ2)
      call myBget(p1Pp22,MQ2,MQ2)

      color=-1d0/(2d0*Nc)
      Ap=Ap+color*(GS**3*(2*MGl*MQ*mixing*(xc1 - xc2) + (xb0 + 2*xb1 - MGl2*xc1 - MQ2*xc1 + p1p1*xc1 + MGl2*xc2 + MQ2*xc2 - p2p2*xc2)*Krondelta(sf1,sf2)))/(8d0*Pi**2)
      Am=Am-color*(GS**3*(-2*MGl*MQ*mixing*(xc0 + xc1 + xc2)+ (xb0 + 2*MGl2*xc0 + MGl2*xc1 + MQ2*xc1 + p1p1*xc1 + MGl2*xc2 + MQ2*xc2 + p2p2*xc2)*Krondelta(sf1,sf2)))/(8d0*Pi**2)

! Quark exchange (8) (Quark-Gluino-Gluino)
      call myCget(p1p1,p1Pp22,p2p2,MQ2,MGl2,MGl2)
      call myBget(p1Pp22,MGl2,MGl2)

      color=Nc/2d0
      Ap=Ap+color*(GS**3*(2*MGl*MQ*mixing*(xc1 - xc2) + (xb0 + 2*xb1- MGl2*xc1 - MQ2*xc1 + p1p1*xc1 + MGl2*xc2 + MQ2*xc2 - p2p2*xc2)*Krondelta(sf1,sf2)))/(8d0*Pi**2)
      Am=Am-color*(GS**3*(-2*MGl*MQ*mixing*(xc0 + xc1 + xc2) + (xb0 + 2*MQ2*xc0 + MGl2*xc1 + MQ2*xc1 + p1p1*xc1 + MGl2*xc2 + MQ2*xc2 + p2p2*xc2)*Krondelta(sf1,sf2)))/(8d0*Pi**2)

      endif

      end



      

! The function "stsT2xx_PiGluonLoopT" computes the transverse part of the gluon propagator correction.
      double complex function stsT2xx_PiGluonLoopT(myP2)
      implicit none
#include "stsT2xx_Kinematics.h"
#include "stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer i,k,m
      double precision ML02,myP2
      
      stsT2xx_PiGluonLoopT=0d0

! ===== Quark loop =====  I
      do k=3,4
      do m=1,3

      ML02 = Mf2(k,m)

      if (m.le.2) then
      ML02=0d0
      endif

      call myBget(myP2,ML02,ML02)
      call myAget(ML02)

      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT-GS**2*TF/(8d0*pi**2)*(2d0*xa0-4d0*xb00-myP2*xb0)        

      enddo
      enddo

! ===== Gluino loop ===== II

      call myBget(myP2,MGl2,MGl2)
      call myAget(MGl2)

      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT-GS**2*CV/(16d0*pi**2)*(2d0*xa0-4d0*xb00-myP2*xb0)        
        
        
! ===== Gluon loop =====  III
      call myBget(myP2,0d0,0d0)

      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT-GS**2*CV/(16d0*pi**2)*(2d0*myP2*xb0+5d0*xb00)

! ===== Squark loop =====     IV

      do k=3,4
      do m=1,3
      do i=1,2

      ML02 = MSf2(i,k,m)

      call myBget(myP2,ML02,ML02)

      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT-GS**2*TF/(4d0*pi**2)*xb00

      enddo
      enddo
      enddo 

! ===== Ghost loop =====  V

      call myBget(myP2,0d0,0d0)

      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT+GS**2*CV/(16d0*pi**2)*xb00

! ----------Gluon-onepoint-Loop vanishes-----------  VI
! --------------Squark-onepoint-Loop-----------------------   VII

      do i=1,2
      do k=3,4
      do m=1,3
      ML02 = MSf2(i,k,m)
      call myAget(ML02)
      stsT2xx_PiGluonLoopT=stsT2xx_PiGluonLoopT+GS**2*TF/(8d0*Pi**2)*xa0
      enddo
      enddo
      enddo   
       
      end