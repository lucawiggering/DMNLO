      subroutine stst2gg_calc_A4(PA_A4,i,type1,gen1,j,type2,gen2)
      implicit none 
!i,type1,gen1 are sfermion index, type and generation of squark
!j,type2,gen2 are sfermion index, type and generation of antisquark
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"
#include "../common/stsT2xx_Kinematics.h"

!     The first index of PA_A4 (PartialAmplitudes) denotes the colour basis vector: 
!     1---> singlet (1), 2---> symmetric octet (8_S), 3---> antisymmetric octet (8_A) 
!     The second index denotes the kinetic factor: 
!     1---> g^{\mu\nu}, 
!     2--> k_1^\mu k_1^\nu, 3--> k_2^\mu k_2^\nu, 4--> p_1^\nu p_1^\mu, 5--> p_2^\nu p_2^\mu 
!     6--> k_1^\mu k_2^\nu, 7--> k_1^\mu p_1^\nu, 8--> k_1^\mu p_2^\nu,
!     9--> k_2^\mu k_1^\nu, 10--> k_2^\mu p_1^\nu, 11--> k_2^\mu p_2^\nu,
!     12--> p_1^\mu k_1^\nu, 13--> p_1^\mu k_2^\nu, 14--> p_1^\mu p_2^\nu,
!     15--> p_2^\mu k_1^\nu, 16--> p_2^\mu k_2^\nu, 17--> p_2^\mu p_1^\nu,


      double complex PA_A4(3,17), temp,tempMixing,tempPA(17),c0
      integer i,j,k,type1,type2,gen1,gen2,n,m,gammaType,gammaGen
      double precision colors(3),M2,M2j,MQ,MQ2,Krondelta
      
      c0=(0d0,0d0)
      PA_A4=0d0
      M2=MSf2(i,type1,gen1)
      M2j=MSf2(j,type2,gen2)
      
      if (gen1.eq.3) then 
      MQ=Mf(type1,gen1)
      MQ2=Mf2(type1,gen1)
      else
      MQ=0d0
      MQ2=0d0
      endif

!     The first condition excludes the annihilation of squarks with different flavours.
      if (type1.eq.type2 .and. gen1.eq.gen2) then
!     The the second condition excludes all corrections which allow for different sfermion numbers. 
!     The SquarkLoop (1), SquarkExchangeRight (16) and Boxes 4-6 allow for i!=j. These are given at the end.
      if (i.eq.j) then

!     Gluon Loop (2) 
      call myBget(sman,0d0,0d0)

      colors=(/1d0,Nc/2d0,0d0/)
      temp=GS**4/(16d0*Pi**2)*(-3*xb0)

      do m=1,2
      PA_A4(m,1)=PA_A4(m,1)+colors(m)*temp
      enddo

!     Squark-Gluon Loop (3) 
      call myBget(tman,0d0,M2)

      colors=(/(Nc**2-2d0)/(2d0*Nc**2),(Nc**2-4d0)/(4d0*Nc),-(Nc**2-4d0)/(4d0*Nc)/)
      temp=-GS**4/(16d0*Pi**2)*xb0

      do m=1,3
      PA_A4(m,1)=PA_A4(m,1)+colors(m)*temp
      enddo

!     Squark-Gluon Loop u (4) 
      call myBget(uman,0d0,M2)

      colors=(/(Nc**2-2d0)/(2d0*Nc**2),(Nc**2-4d0)/(4d0*Nc),(Nc**2-4d0)/(4d0*Nc)/)
      temp=-GS**4/(16d0*Pi**2)*xb0

      do m=1,3
      PA_A4(m,1)=PA_A4(m,1)+colors(m)*temp
      enddo
      
!     GluonExchangeLeft (5) 
      call myCget(M2,sman,M2,0d0,M2,M2)
      call myBget(sman,M2,M2)

      temp=(GS**4)/(16d0*Pi**2)*(xb0+(4*M2-2*sman)*xc0+2*(4*M2-sman)*xc1)
      PA_A4(1,1)=PA_A4(1,1)+CF/Nc*temp
      PA_A4(2,1)=PA_A4(2,1)-1d0/(2d0*Nc)*temp

!     GluonExchangeAbove (6) 
      call myCget(M2,tman,0d0,0d0,M2,0d0) ! double poles
      call myBget(tman,M2,0d0)

      colors=(/0d0,-Nc/4d0,Nc/4d0/)
      tempPA=(/-(GS**4*(xb0-2*tman*xc0-xc00-tman*xc1-tman*xc2+M2*(2*xc0+3*xc1+xc2)))/(16d0*Pi**2),
     &(GS**4*(2*xc2+xc22))/(16d0*Pi**2),c0,
     &(GS**4*(2*xc1+xc11))/(16d0*Pi**2),c0,c0,
     &(GS**4*(2*xc0+xc1+xc12+4*xc2))/(16d0*Pi**2),c0,c0,c0,c0,
     &(GS**4*(2*xc0+xc1+xc12-2*xc2))/(16d0*Pi**2),c0,c0,c0,c0,c0/)

      do n=1,12
      do m=2,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo
     
!     GluonExchangeAboveU (7)
      call myCget(M2,uman,0d0,0d0,M2,0d0) !double poles
      call myBget(uman,M2,0d0)

      colors=(/0d0,-Nc/4d0,-Nc/4d0/)
      tempPA=(/-(GS**4*(xb0-2*uman*xc0-xc00-uman*xc1-uman*xc2+M2*(2*xc0+3*xc1+xc2)))/(16d0*Pi**2),c0,
     &(GS**4*(2*xc2+xc22))/(16d0*Pi**2),
     &(GS**4*(2*xc1+xc11))/(16d0*Pi**2),c0,c0,c0,c0,c0,
     &(GS**4*(2*xc0+xc1+xc12-2*xc2))/(16d0*Pi**2),c0,c0,
     &(GS**4*(2*xc0+xc1+xc12+4*xc2))/(16d0*Pi**2),c0,c0,c0,c0/)

      do n=1,13
      do m=2,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     GluonExchangeBelow (8) 
      call myCget(M2,tman,0d0,0d0,M2,0d0)!double poles
      call myBget(tman,M2,0d0)

      colors=(/0d0,Nc/4d0,-Nc/4d0/)
      tempPA=(/(GS**4*(xb0-2*tman*xc0-xc00-tman*xc1-tman*xc2+M2*(2*xc0+3*xc1+xc2)))/(16d0*Pi**2),c0,
     &-(GS**4*(2*xc2+xc22))/(16d0*Pi**2),c0,
     &-(GS**4*(2*xc1+xc11))/(16d0*Pi**2),c0,c0,c0,c0,c0,
     &-(GS**4*(2*xc0+xc1+xc12-2*xc2))/(16d0*Pi**2),c0,c0,c0,c0,
     &-(GS**4*(2*xc0+xc1+xc12+4*xc2))/(16d0*Pi**2),c0/)

      do n=1,17
      do m=2,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo
   
!     GluonExchangeBelowU (9) 
      call myCget(M2,uman,0d0,0d0,M2,0d0) !double poles
      call myBget(uman,M2,0d0)

      colors=(/0d0,Nc/4d0,Nc/4d0/)
      tempPA=(/(GS**4*(xb0-2*uman*xc0-xc00-uman*xc1-uman*xc2+M2*(2*xc0+3*xc1+xc2)))/(16d0*Pi**2),
     &-(GS**4*(2*xc2+xc22))/(16d0*Pi**2),c0,c0,
     &-(GS**4*(2*xc1+xc11))/(16d0*Pi**2),c0,c0,
     &-(GS**4*(2*xc0+xc1+xc12+4*xc2))/(16d0*Pi**2),c0,c0,c0,c0,c0,c0,
     &-(GS**4*(2*xc0+xc1+xc12-2*xc2))/(16d0*Pi**2),c0,c0/)

      do n=1,17
      do m=2,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     GluonExchangeRight (10) 
      call myCget(0d0,sman,0d0,0d0,0d0,0d0) ! double poles/ double checked
      call myBget(sman,0d0,0d0)

      colors=(/1d0,Nc/2d0,0d0/)
      tempPA=(/(GS**4)/(16d0*Pi**2)*(2*xb0-5d0/2d0*sman*xc0+10*xc00-sman*xc1),
     &(GS**4)/(16d0*Pi**2)*(5*xc1+10*xc11),(GS**4)/(16d0*Pi**2)*(5*xc1+10*xc11),c0,c0,(-GS**4)/(8d0*Pi**2)*(xc0+4*xc1+5*xc12),
     &c0,c0,(GS**4)/(8d0*Pi**2)*(2*xc0-xc1-5*xc12),c0,c0,c0,c0,c0,c0,c0,c0/)
      
      do n=1,9
      do m=1,2
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     SquarkExchangeLeft (11) 
      call myCget(M2,sman,M2,M2,0d0,0d0)
      call myBget(sman,0d0,0d0)

      colors=(/1d0,Nc/2d0,0d0/)
      tempPA=(/(GS**4)/(32d0*Pi**2)*(2*xb0+(4*M2-sman)*xc0-2*xc00+2*(sman-4*M2)*xc1),c0,c0,
     &(GS**4)/(16d0*Pi**2)*(xc1-xc11),(GS**4)/(16d0*Pi**2)*(xc1-xc11), c0,c0,c0,c0,c0,c0,c0,
     &c0,(GS**4)/(32d0*Pi**2)*(xc0-2*xc1+2*xc12),c0,c0,(GS**4)/(32d0*Pi**2)*(xc0-2*xc1+2*xc12)/)

      do n=1,17
      do m=1,2
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

      temp=3*Nc*GS**4/(64d0*Pi**2)*(xc0-2*xc1)
      PA_A4(3,14)=PA_A4(3,14)+temp
      PA_A4(3,17)=PA_A4(3,17)-temp

!     SquarkExchangeAbove (12) 
      call myCget(0d0,M2,tman,M2,M2,0d0)

      colors=(/(Nc**2-2d0)/(4d0*Nc**2),-1d0/(2d0*Nc),1d0/(2d0*Nc)/)
      tempPA=(/(GS**4)/(8d0*Pi**2)*xc00,(GS**4)/(16d0*Pi**2)*(xc1+2*xc11),c0,(-GS**4)/(8d0*Pi**2)*(xc2-xc22),c0,c0,
     &(-GS**4)/(16d0*Pi**2)*(xc0+2*xc1-2*xc12-xc2),c0,c0,c0,c0,(GS**4)/(8d0*Pi**2)*xc12,c0,c0,c0,c0,c0/)

      do n=1,12
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     SquarkExchangeAboveU (13) 
      call myCget(0d0,M2,uman,M2,M2,0d0) ! no double poles

      colors=(/(Nc**2-2d0)/(4d0*Nc**2),-1d0/(2d0*Nc),-1d0/(2d0*Nc)/)
      tempPA=(/(GS**4)/(8d0*Pi**2)*xc00,c0,(GS**4)/(16d0*Pi**2)*(xc1+2*xc11),(-GS**4)/(8d0*Pi**2)*(xc2-xc22),c0,c0,
     &c0,c0,c0,(GS**4)/(8d0*Pi**2)*xc12,c0,c0,(-GS**4)/(16d0*Pi**2)*(xc0+2*xc1-2*xc12-xc2),c0,c0,c0,c0/)

      do n=1,13
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     SquarkExchangeBelow (14) 
      call myCget(0d0,tman,M2,M2,M2,0d0) ! no double poles

      colors=(/(Nc**2-2d0)/(4d0*Nc**2),-1d0/(2d0*Nc),1d0/(2d0*Nc)/)
      tempPA=(/(GS**4)/(8d0*Pi**2)*xc00,c0,(GS**4)/(16d0*Pi**2)*(xc1+2*xc11),c0,(-GS**4)/(8d0*Pi**2)*(xc2-xc22),c0,c0,
     &c0,c0,c0,(GS**4)/(8d0*Pi**2)*xc12,c0,c0,c0,c0,(-GS**4)/(16d0*Pi**2)*(xc0+2*xc1-2*xc12-xc2),c0/)

      do n=1,16
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     SquarkExchangeBelowU (15) ! no double poles
      call myCget(0d0,uman,M2,M2,M2,0d0)

      colors=(/(Nc**2-2d0)/(4d0*Nc**2),-1d0/(2d0*Nc),-1d0/(2d0*Nc)/)
      tempPA=(/(GS**4)/(8d0*Pi**2)*xc00,(GS**4)/(16d0*Pi**2)*(xc1+2*xc11),c0,c0,(-GS**4)/(8d0*Pi**2)*(xc2-xc22),
     &c0,c0,(-GS**4)/(16d0*Pi**2)*(xc0+2*xc1-2*xc12-xc2),c0,c0,c0,c0,c0,c0,(GS**4)/(8d0*Pi**2)*xc12,c0,c0/)

      do n=1,15
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     Box 1 (17) AGREES WITH FORMCALC
      call LTDget(M2,0d0,0d0,M2,tman,sman,0d0,M2,M2,M2) 
      call myCget(0d0,0d0,sman,M2,M2,M2)

      colors=(/CF/(2d0*Nc),-1d0/(4d0*Nc),1d0/(4d0*Nc)/)
      tempPA=(/-(GS**4*(xc00-2*sman*xd00-sman*xd001-sman*xd002+tman*xd002-uman*xd002-sman*xd003+4*M2*(xd00+xd001+xd002+xd003)))/(4d0*Pi**2),c0,
     &-(GS**4*(xc0+4*xc1+4*xc11+4*M2*xd0-2*sman*xd0+4*M2*xd1-sman*xd1+16*M2*xd12-4*sman*xd12+16*M2*xd122-4*sman*xd122+20*M2*xd2-9*sman*xd2+tman*xd2-uman*xd2+32*M2*xd22-12*sman*xd22+4*tman*xd22-4*uman*xd22+16*M2*xd222-4*sman*xd222+4*tman*xd222-4*uman*xd222+16*M2*xd223-4*sman*xd223+16*M2*xd23-4*sman*xd23+4*M2*xd3-sman*xd3))/(16d0*Pi**2),
     &-(GS**4*(xc2+2*xc22+2*xd00+8*xd003+4*M2*xd13-sman*xd13+8*M2*xd133-2*sman*xd133+4*M2*xd23-sman*xd23+tman*xd23-uman*xd23+8*M2*xd233-2*sman*xd233+2*tman*xd233-2*uman*xd233+4*M2*xd3-2*sman*xd3+12*M2*xd33-5*sman*xd33+8*M2*xd333-2*sman*xd333))/(8d0*Pi**2),
     &-(GS**4*(xc2+2*xc22+4*M2*xd0-2*sman*xd0+6*xd00+8*xd001+8*xd002+16*M2*xd1-7*sman*xd1+20*M2*xd11-7*sman*xd11+8*M2*xd111-2*sman*xd111+24*M2*xd112-6*sman*xd112+2*tman*xd112-2*uman*xd112+8*M2*xd113-2*sman*xd113+40*M2*xd12-14*sman*xd12+3*tman*xd12-3*uman*xd12+24*M2*xd122-6*sman*xd122+4*tman*xd122-4*uman*xd122+16*M2*xd123-4*sman*xd123+12*M2*xd13-3*sman*xd13+16*M2*xd2-7*sman*xd2+tman*xd2-uman*xd2+20*M2*xd22-7*sman*xd22+3*tman*xd22-3*uman*xd22+8*M2*xd222-2*sman*xd222+2*tman*xd222-2*uman*xd222+8*M2*xd223-2*sman*xd223+12*M2*xd23-3*sman*xd23+4*M2*xd3-sman*xd3))/(8d0*Pi**2),c0,c0,c0,c0,
     &-(GS**4*(2*xc12+xc2+2*xd00+4*xd002+8*M2*xd123-2*sman*xd123+4*M2*xd13-sman*xd13+8*M2*xd223-2*sman*xd223+2*tman*xd223-2*uman*xd223+12*M2*xd23-5*sman*xd23+tman*xd23-uman*xd23+8*M2*xd233-2*sman*xd233+4*M2*xd3-2*sman*xd3+4*M2*xd33-sman*xd33))/(8d0*Pi**2),
     &-(GS**4*(2*xc12+xc2-4*M2*xd0+2*sman*xd0-2*xd00-4*xd002-8*M2*xd1+3*sman*xd1-4*M2*xd11+sman*xd11-8*M2*xd112+2*sman*xd112-24*M2*xd12+8*sman*xd12-tman*xd12+uman*xd12-16*M2*xd122+4*sman*xd122-2*tman*xd122+2*uman*xd122-8*M2*xd123+2*sman*xd123-4*M2*xd13+sman*xd13-16*M2*xd2+7*sman*xd2-tman*xd2+uman*xd2-20*M2*xd22+7*sman*xd22-3*tman*xd22+3*uman*xd22-8*M2*xd222+2*sman*xd222-2*tman*xd222+2*uman*xd222-8*M2*xd223+2*sman*xd223-12*M2*xd23+3*sman*xd23-4*M2*xd3+sman*xd3))/(8d0*Pi**2),c0,
     &-(GS**4*(xc0+2*xc1+4*xc12+2*xc2+4*M2*xd0-2*sman*xd0+4*xd00+8*xd002+4*M2*xd1-sman*xd1+8*M2*xd12-2*sman*xd12+16*M2*xd123-4*sman*xd123+8*M2*xd13-2*sman*xd13+12*M2*xd2-5*sman*xd2+tman*xd2-uman*xd2+8*M2*xd22-2*sman*xd22+2*tman*xd22-2*uman*xd22+16*M2*xd223-4*sman*xd223+4*tman*xd223-4*uman*xd223+32*M2*xd23-12*sman*xd23+2*tman*xd23-2*uman*xd23+16*M2*xd233-4*sman*xd233+12*M2*xd3-5*sman*xd3+8*M2*xd33-2*sman*xd33))/(16d0*Pi**2),
     &-(GS**4*(xc2+2*xc22-4*M2*xd0+2*sman*xd0-6*xd00-4*xd001-4*xd002-4*xd003-8*M2*xd1+3*sman*xd1-4*M2*xd11+sman*xd11-8*M2*xd113+2*sman*xd113-8*M2*xd12+2*sman*xd12-tman*xd12+uman*xd12-16*M2*xd123+4*sman*xd123-2*tman*xd123+2*uman*xd123-20*M2*xd13+7*sman*xd13-8*M2*xd133+2*sman*xd133-8*M2*xd2+3*sman*xd2-tman*xd2+uman*xd2-4*M2*xd22+sman*xd22-tman*xd22+uman*xd22-8*M2*xd223+2*sman*xd223-2*tman*xd223+2*uman*xd223-20*M2*xd23+7*sman*xd23-2*tman*xd23+2*uman*xd23-8*M2*xd233+2*sman*xd233-12*M2*xd3+5*sman*xd3-8*M2*xd33+2*sman*xd33))/(8d0*Pi**2),c0,
     &-(GS**4*(xc0+2*xc1+4*xc12+2*xc2-4*M2*xd0+2*sman*xd0-4*xd00-8*xd002-12*M2*xd1+5*sman*xd1-8*M2*xd11+2*sman*xd11-16*M2*xd112+4*sman*xd112-40*M2*xd12+14*sman*xd12-2*tman*xd12+2*uman*xd12-32*M2*xd122+8*sman*xd122-4*tman*xd122+4*uman*xd122-16*M2*xd123+4*sman*xd123-8*M2*xd13+2*sman*xd13-20*M2*xd2+9*sman*xd2-tman*xd2+uman*xd2-32*M2*xd22+12*sman*xd22-4*tman*xd22+4*uman*xd22-16*M2*xd222+4*sman*xd222-4*tman*xd222+4*uman*xd222-16*M2*xd223+4*sman*xd223-16*M2*xd23+4*sman*xd23-4*M2*xd3+sman*xd3))/(16d0*Pi**2),
     &-(GS**4*(xc2+2*xc22-2*xd00-4*xd001-4*xd002-4*xd003-8*M2*xd113+2*sman*xd113-16*M2*xd123+4*sman*xd123-2*tman*xd123+2*uman*xd123-12*M2*xd13+5*sman*xd13-8*M2*xd133+2*sman*xd133-8*M2*xd223+2*sman*xd223-2*tman*xd223+2*uman*xd223-12*M2*xd23+5*sman*xd23-tman*xd23+uman*xd23-8*M2*xd233+2*sman*xd233-4*M2*xd3+2*sman*xd3-4*M2*xd33+sman*xd33))/(8d0*Pi**2)/)
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     Box 1 u (18) AGREES WITH FORMCALC
      call LTDget(M2,0d0,0d0,M2,uman,sman,0d0,M2,M2,M2) 
      call myCget(0d0,0d0,sman,M2,M2,M2)

      colors=(/CF/(2d0*Nc),-1d0/(4d0*Nc),-1d0/(4d0*Nc)/)
      tempPA=(/-(GS**4*(xc00-2*sman*xd00-sman*xd001-sman*xd002-tman*xd002+uman*xd002-sman*xd003+4*M2*(xd00+xd001+xd002+xd003)))/(4d0*Pi**2),
     &-(GS**4*(xc0+4*xc1+4*xc11+4*M2*xd0-2*sman*xd0+4*M2*xd1-sman*xd1+16*M2*xd12-4*sman*xd12+16*M2*xd122-4*sman*xd122+20*M2*xd2-9*sman*xd2-tman*xd2+uman*xd2+32*M2*xd22-12*sman*xd22-4*tman*xd22+4*uman*xd22+16*M2*xd222-4*sman*xd222-4*tman*xd222+4*uman*xd222+16*M2*xd223-4*sman*xd223+16*M2*xd23-4*sman*xd23+4*M2*xd3-sman*xd3))/(16d0*Pi**2),c0,
     &-(GS**4*(xc2+2*xc22+2*xd00+8*xd003+4*M2*xd13-sman*xd13+8*M2*xd133-2*sman*xd133+4*M2*xd23-sman*xd23-tman*xd23+uman*xd23+8*M2*xd233-2*sman*xd233-2*tman*xd233+2*uman*xd233+4*M2*xd3-2*sman*xd3+12*M2*xd33-5*sman*xd33+8*M2*xd333-2*sman*xd333))/(8d0*Pi**2),
     &-(GS**4*(xc2+2*xc22+4*M2*xd0-2*sman*xd0+6*xd00+8*xd001+8*xd002+16*M2*xd1-7*sman*xd1+20*M2*xd11-7*sman*xd11+8*M2*xd111-2*sman*xd111+24*M2*xd112-6*sman*xd112-2*tman*xd112+2*uman*xd112+8*M2*xd113-2*sman*xd113+40*M2*xd12-14*sman*xd12-3*tman*xd12+3*uman*xd12+24*M2*xd122-6*sman*xd122-4*tman*xd122+4*uman*xd122+16*M2*xd123-4*sman*xd123+12*M2*xd13-3*sman*xd13+16*M2*xd2-7*sman*xd2-tman*xd2+uman*xd2+20*M2*xd22-7*sman*xd22-3*tman*xd22+3*uman*xd22+8*M2*xd222-2*sman*xd222-2*tman*xd222+2*uman*xd222+8*M2*xd223-2*sman*xd223+12*M2*xd23-3*sman*xd23+4*M2*xd3-sman*xd3))/(8d0*Pi**2),c0,
     &-(GS**4*(xc0+2*xc1+4*xc12+2*xc2+4*M2*xd0-2*sman*xd0+4*xd00+8*xd002+4*M2*xd1-sman*xd1+8*M2*xd12-2*sman*xd12+16*M2*xd123-4*sman*xd123+8*M2*xd13-2*sman*xd13+12*M2*xd2-5*sman*xd2-tman*xd2+uman*xd2+8*M2*xd22-2*sman*xd22-2*tman*xd22+2*uman*xd22+16*M2*xd223-4*sman*xd223-4*tman*xd223+4*uman*xd223+32*M2*xd23-12*sman*xd23-2*tman*xd23+2*uman*xd23+16*M2*xd233-4*sman*xd233+12*M2*xd3-5*sman*xd3+8*M2*xd33-2*sman*xd33))/(16d0*Pi**2),
     &-(GS**4*(xc0+2*xc1+4*xc12+2*xc2-4*M2*xd0+2*sman*xd0-4*xd00-8*xd002-12*M2*xd1+5*sman*xd1-8*M2*xd11+2*sman*xd11-16*M2*xd112+4*sman*xd112-40*M2*xd12+14*sman*xd12+2*tman*xd12-2*uman*xd12-32*M2*xd122+8*sman*xd122+4*tman*xd122-4*uman*xd122-16*M2*xd123+4*sman*xd123-8*M2*xd13+2*sman*xd13-20*M2*xd2+9*sman*xd2+tman*xd2-uman*xd2-32*M2*xd22+12*sman*xd22+4*tman*xd22-4*uman*xd22-16*M2*xd222+4*sman*xd222+4*tman*xd222-4*uman*xd222-16*M2*xd223+4*sman*xd223-16*M2*xd23+4*sman*xd23-4*M2*xd3+sman*xd3))/(16d0*Pi**2),c0,c0,c0,
     &-(GS**4*(2*xc12+xc2+2*xd00+4*xd002+8*M2*xd123-2*sman*xd123+4*M2*xd13-sman*xd13+8*M2*xd223-2*sman*xd223-2*tman*xd223+2*uman*xd223+12*M2*xd23-5*sman*xd23-tman*xd23+uman*xd23+8*M2*xd233-2*sman*xd233+4*M2*xd3-2*sman*xd3+4*M2*xd33-sman*xd33))/(8d0*Pi**2),c0,
     &-(GS**4*(xc2+2*xc22-2*xd00-4*xd001-4*xd002-4*xd003-8*M2*xd113+2*sman*xd113-16*M2*xd123+4*sman*xd123+2*tman*xd123-2*uman*xd123-12*M2*xd13+5*sman*xd13-8*M2*xd133+2*sman*xd133-8*M2*xd223+2*sman*xd223+2*tman*xd223-2*uman*xd223-12*M2*xd23+5*sman*xd23+tman*xd23-uman*xd23-8*M2*xd233+2*sman*xd233-4*M2*xd3+2*sman*xd3-4*M2*xd33+sman*xd33))/(8d0*Pi**2),
     &-(GS**4*(2*xc12+xc2-4*M2*xd0+2*sman*xd0-2*xd00-4*xd002-8*M2*xd1+3*sman*xd1-4*M2*xd11+sman*xd11-8*M2*xd112+2*sman*xd112-24*M2*xd12+8*sman*xd12+tman*xd12-uman*xd12-16*M2*xd122+4*sman*xd122+2*tman*xd122-2*uman*xd122-8*M2*xd123+2*sman*xd123-4*M2*xd13+sman*xd13-16*M2*xd2+7*sman*xd2+tman*xd2-uman*xd2-20*M2*xd22+7*sman*xd22+3*tman*xd22-3*uman*xd22-8*M2*xd222+2*sman*xd222+2*tman*xd222-2*uman*xd222-8*M2*xd223+2*sman*xd223-12*M2*xd23+3*sman*xd23-4*M2*xd3+sman*xd3))/(8d0*Pi**2),c0,
     &-(GS**4*(xc2+2*xc22-4*M2*xd0+2*sman*xd0-6*xd00-4*xd001-4*xd002-4*xd003-8*M2*xd1+3*sman*xd1-4*M2*xd11+sman*xd11-8*M2*xd113+2*sman*xd113-8*M2*xd12+2*sman*xd12+tman*xd12-uman*xd12-16*M2*xd123+4*sman*xd123+2*tman*xd123-2*uman*xd123-20*M2*xd13+7*sman*xd13-8*M2*xd133+2*sman*xd133-8*M2*xd2+3*sman*xd2+tman*xd2-uman*xd2-4*M2*xd22+sman*xd22+tman*xd22-uman*xd22-8*M2*xd223+2*sman*xd223+2*tman*xd223-2*uman*xd223-20*M2*xd23+7*sman*xd23+2*tman*xd23-2*uman*xd23-8*M2*xd233+2*sman*xd233-12*M2*xd3+5*sman*xd3-8*M2*xd33+2*sman*xd33))/(8d0*Pi**2)/)
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     Box 2 (19) A23 (11) AND A43 (17) DO NOT AGREE WITH FORMCALC
      call LTDget(M2,0d0,0d0,M2,tman,sman,M2,0d0,0d0,0d0)  
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
      call myBget(0d0,0d0,0d0)

      colors=(/1d0/2d0,Nc/4d0,-Nc/4d0/)
      tempPA=(/-(GS**4*(xb0-3*tman*xc0+uman*xc0-xc00+uman*xc1-sman*xc2-2*tman*xc2-2*uman*xc2-sman*tman*xd0-tman*uman*xd0-7*sman*xd00+2*tman*xd00+3*uman*xd00+sman*xd001-tman*xd001+sman*xd002-3*tman*xd002+uman*xd002+2*sman*xd003+uman*xd003-tman**2*xd1+tman*uman*xd1+sman*tman*xd11+tman**2*xd11+2*sman*tman*xd12+3*tman**2*xd12+tman*uman*xd12-sman*tman*xd13-sman*uman*xd13-2*tman*uman*xd13-2*tman**2*xd2+sman*tman*xd22+2*tman**2*xd22+tman*uman*xd22-sman*tman*xd23-sman*uman*xd23-3*tman*uman*xd23-uman**2*xd23+sman*tman*xd3-sman*uman*xd3+tman*uman*xd3-uman**2*xd3+sman*uman*xd33+uman**2*xd33+M2**2*(-xd0+xd11+4*xd12-2*xd13-2*xd2+3*xd22-4*xd23+xd33)+M2*(2*xc0-xc1+4*xc2+sman*xd0+tman*xd0+uman*xd0+5*xd00-5*xd001-4*xd002-7*xd003+tman*xd1
     &-uman*xd1-sman*xd11-2*tman*xd11-2*sman*xd12-7*tman*xd12-uman*xd12+2*sman*xd13+2*tman*xd13+2*uman*xd13+4*tman*xd2-sman*xd22-5*tman*xd22-uman*xd22+2*sman*xd23+3*tman*xd23+5*uman*xd23-tman*xd3+uman*xd3-sman*xd33-2*uman*xd33)))/(16d0*Pi**2),c0,
     &-(GS**4*(2*xc0-4*xc1-2*xc11+4*M2*xd0-sman*xd0-12*xd00-8*xd002-4*M2*xd1+sman*xd1-10*M2*xd12+sman*xd12-6*tman*xd12-10*M2*xd122+2*sman*xd122-2*tman*xd122+12*M2*xd2-10*sman*xd2+2*tman*xd2-2*uman*xd2+14*M2*xd22-13*sman*xd22-13*tman*xd22+3*uman*xd22-8*M2*xd222+2*sman*xd222-6*tman*xd222+2*uman*xd222-14*M2*xd223+4*sman*xd223+2*uman*xd223-22*M2*xd23+7*sman*xd23+6*uman*xd23-4*M2*xd3+sman*xd3))/(32d0*Pi**2),
     &(GS**4*(2*xc0-8*xc2+2*xc22+2*M2*xd0-2*tman*xd0-2*xd00+8*xd003-2*M2*xd1+2*tman*xd1+17*M2*xd13+sman*xd13-5*tman*xd13+10*M2*xd133-2*sman*xd133+2*tman*xd133-2*M2*xd2+2*tman*xd2+9*M2*xd23+sman*xd23+2*tman*xd23+uman*xd23+8*M2*xd233-2*sman*xd233+6*tman*xd233-2*uman*xd233-5*M2*xd3-3*sman*xd3+2*tman*xd3-uman*xd3-11*M2*xd33+7*sman*xd33+3*uman*xd33+14*M2*xd333-4*sman*xd333-2*uman*xd333))/(32d0*Pi**2),
     &-(GS**4*(2*xc0-8*xc2-2*xc22+7*M2*xd0+sman*xd0+uman*xd0-2*xd00-16*xd001-16*xd002+12*M2*xd1-2*sman*xd1+3*tman*xd1-3*uman*xd1-9*M2*xd11-sman*xd11-tman*xd11+2*uman*xd11-10*M2*xd111+2*sman*xd111-2*tman*xd111-28*M2*xd112+6*sman*xd112-10*tman*xd112+2*uman*xd112-14*M2*xd113+4*sman*xd113+2*uman*xd113-30*M2*xd12-2*sman*xd12-3*tman*xd12+17*uman*xd12-26*M2*xd122+6*sman*xd122-14*tman*xd122+4*uman*xd122-28*M2*xd123+8*sman*xd123+4*uman*xd123-29*M2*xd13+sman*xd13+uman*xd13+6*M2*xd2-2*sman*xd2+4*tman*xd2+2*uman*xd2-21*M2*xd22-sman*xd22-2*tman*xd22+15*uman*xd22-8*M2*xd222+2*sman*xd222-6*tman*xd222+2*uman*xd222-14*M2*xd223+4*sman*xd223+2*uman*xd223-29*M2*xd23+sman*xd23+uman*xd23-7*M2*xd3-sman*xd3-uman*xd3))/(32d0*Pi**2),c0,c0,c0,c0,
     &-(GS**4*(2*xc0+8*xc1-2*xc12-2*xc2+2*M2*xd0-2*tman*xd0+18*xd00-4*xd002-4*xd003-2*M2*xd1+2*tman*xd1-8*M2*xd12+8*tman*xd12-10*M2*xd123+2*sman*xd123-2*tman*xd123-5*M2*xd13-sman*xd13-3*tman*xd13+6*M2*xd2-6*tman*xd2-8*M2*xd22+8*tman*xd22-8*M2*xd223+2*sman*xd223-6*tman*xd223+2*uman*xd223+27*M2*xd23-15*sman*xd23-8*tman*xd23-7*uman*xd23-14*M2*xd233+4*sman*xd233+2*uman*xd233+9*M2*xd3-5*sman*xd3+2*tman*xd3-3*uman*xd3-11*M2*xd33+5*sman*xd33+3*uman*xd33))/(32d0*Pi**2),
     &(GS**4*(10*xc0+2*xc12+2*xc2+7*M2*xd0+sman*xd0+uman*xd0-22*xd00-4*xd001-12*xd002-2*M2*xd1+3*tman*xd1-uman*xd1-5*M2*xd11-sman*xd11-3*tman*xd11-10*M2*xd112+2*sman*xd112-2*tman*xd112-10*sman*xd12-17*tman*xd12+uman*xd12-18*M2*xd122+4*sman*xd122-8*tman*xd122+2*uman*xd122-14*M2*xd123+4*sman*xd123+2*uman*xd123-11*M2*xd13+5*sman*xd13+3*uman*xd13+14*M2*xd2-10*sman*xd2+4*tman*xd2-6*uman*xd2+3*M2*xd22-9*sman*xd22-18*tman*xd22+7*uman*xd22-8*M2*xd222+2*sman*xd222-6*tman*xd222+2*uman*xd222-14*M2*xd223+4*sman*xd223+2*uman*xd223-37*M2*xd23+9*sman*xd23+9*uman*xd23-7*M2*xd3-sman*xd3-uman*xd3))/(32d0*Pi**2),c0,
     &-(GS**4*(6*xc0-2*xc12-2*xc2+5*M2*xd0+sman*xd0-tman*xd0-10*xd00-4*xd002-4*xd003-5*M2*xd1-sman*xd1+tman*xd1-9*M2*xd12-sman*xd12-3*tman*xd12-10*M2*xd123+2*sman*xd123-2*tman*xd123-5*M2*xd13+2*sman*xd13-3*tman*xd13-2*M2*xd2+2*sman*xd2+4*tman*xd2-2*uman*xd2-M2*xd22-sman*xd22-10*tman*xd22-uman*xd22-8*M2*xd223+2*sman*xd223-6*tman*xd223+2*uman*xd223-2*M2*xd23-5*sman*xd23-5*tman*xd23+7*uman*xd23-14*M2*xd233+4*sman*xd233+2*uman*xd233+6*M2*xd3-3*sman*xd3+tman*xd3-3*uman*xd3-11*M2*xd33+2*sman*xd33+3*uman*xd33))/(32d0*Pi**2),
     &(GS**4*(14*xc0+2*xc22+11*M2*xd0+sman*xd0+uman*xd0-26*xd00-4*xd001-4*xd002-8*xd003-2*M2*xd1+3*tman*xd1-uman*xd1-9*M2*xd11-sman*xd11-3*tman*xd11-10*M2*xd113+2*sman*xd113-2*tman*xd113-10*M2*xd12-2*sman*xd12-13*tman*xd12-uman*xd12-18*M2*xd123+4*sman*xd123-8*tman*xd123+2*uman*xd123-11*M2*xd13-sman*xd13-6*tman*xd13+5*uman*xd13-14*M2*xd133+4*sman*xd133+2*uman*xd133+6*M2*xd2-2*tman*xd2-4*uman*xd2-M2*xd22-sman*xd22-10*tman*xd22-uman*xd22-8*M2*xd223+2*sman*xd223-6*tman*xd223+2*uman*xd223-13*M2*xd23-sman*xd23-10*tman*xd23+11*uman*xd23-14*M2*xd233+4*sman*xd233+2*uman*xd233+15*M2*xd3-5*sman*xd3-7*uman*xd3-26*M2*xd33+4*sman*xd33+6*uman*xd33))/(32d0*Pi**2),c0,
     &(GS**4*(6*xc0+8*xc1+2*xc12+2*xc2+2*M2*xd0+sman*xd0+2*uman*xd0+6*xd00-4*xd001-12*xd002+3*M2*xd1-3*sman*xd1+3*tman*xd1-2*uman*xd1-5*M2*xd11+2*sman*xd11-3*tman*xd11-10*M2*xd112+2*sman*xd112-2*tman*xd112-5*M2*xd12-3*sman*xd12-3*tman*xd12+4*uman*xd12-18*M2*xd122+4*sman*xd122-8*tman*xd122+2*uman*xd122-14*M2*xd123+4*sman*xd123+2*uman*xd123-11*M2*xd13+2*sman*xd13+3*uman*xd13-2*M2*xd2+6*tman*xd2+8*uman*xd2-10*M2*xd22-5*sman*xd22+3*tman*xd22+11*uman*xd22-8*M2*xd222+2*sman*xd222-6*tman*xd222+2*uman*xd222-14*M2*xd223+4*sman*xd223+2*uman*xd223-14*M2*xd23-sman*xd23-2*uman*xd23-2*M2*xd3-sman*xd3-2*uman*xd3))/(32d0*Pi**2),
     &(GS**4*(-6*xc0+2*xc22+2*M2*xd0-2*tman*xd0-14*xd00-4*xd001-4*xd002-8*xd003+6*M2*xd1-6*tman*xd1-8*M2*xd11+8*tman*xd11-10*M2*xd113+2*sman*xd113-2*tman*xd113-16*M2*xd12+16*tman*xd12-18*M2*xd123+4*sman*xd123-8*tman*xd123+2*uman*xd123+13*M2*xd13-7*sman*xd13+5*tman*xd13-6*uman*xd13-14*M2*xd133+4*sman*xd133+2*uman*xd133+6*M2*xd2-6*tman*xd2-8*M2*xd22+8*tman*xd22-8*M2*xd223+2*sman*xd223-6*tman*xd223+2*uman*xd223+3*M2*xd23-7*sman*xd23+8*tman*xd23+uman*xd23-14*M2*xd233+4*sman*xd233+2*uman*xd233+M2*xd3+3*sman*xd3+2*tman*xd3+5*uman*xd3-3*M2*xd33-3*sman*xd33-5*uman*xd33))/(32d0*Pi**2)/)      
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     Box 2 u (20) A13 (8) AND A43 (17) DO NOT AGREE WITH FORMCALC
      call LTDget(M2,0d0,0d0,M2,uman,sman,M2,0d0,0d0,0d0)
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
      call myBget(0d0,0d0,0d0)

      colors=(/1d0/2d0,Nc/4d0,Nc/4d0/)
      tempPA=(/-(GS**4*(xb0+tman*xc0-3*uman*xc0-xc00+tman*xc1-sman*xc2-2*tman*xc2-2*uman*xc2-sman*uman*xd0-tman*uman*xd0-7*sman*xd00+3*tman*xd00+2*uman*xd00+sman*xd001-uman*xd001+sman*xd002+tman*xd002-3*uman*xd002+2*sman*xd003+tman*xd003+tman*uman*xd1-uman**2*xd1+sman*uman*xd11+uman**2*xd11+2*sman*uman*xd12+tman*uman*xd12+3*uman**2*xd12-sman*tman*xd13-sman*uman*xd13-2*tman*uman*xd13-2*uman**2*xd2+sman*uman*xd22+tman*uman*xd22+2*uman**2*xd22-sman*tman*xd23-tman**2*xd23-sman*uman*xd23-3*tman*uman*xd23-sman*tman*xd3-tman**2*xd3+sman*uman*xd3+tman*uman*xd3+sman*tman*xd33+tman**2*xd33+M2**2*(-xd0+xd11+4*xd12-2*xd13-2*xd2+3*xd22-4*xd23+xd33)+M2*(2*xc0-xc1+4*xc2+sman*xd0+tman*xd0+uman*xd0+5*xd00-5*xd001-4*xd002
     &-7*xd003-tman*xd1+uman*xd1-sman*xd11-2*uman*xd11-2*sman*xd12-tman*xd12-7*uman*xd12+2*sman*xd13+2*tman*xd13+2*uman*xd13+4*uman*xd2-sman*xd22-tman*xd22-5*uman*xd22+2*sman*xd23+5*tman*xd23+3*uman*xd23+tman*xd3-uman*xd3-sman*xd33-2*tman*xd33)))/(16d0*Pi**2),
     &-(GS**4*(2*xc0-4*xc1-2*xc11+4*M2*xd0-sman*xd0-12*xd00-8*xd002-4*M2*xd1+sman*xd1-10*M2*xd12+sman*xd12-6*uman*xd12-10*M2*xd122+2*sman*xd122-2*uman*xd122+12*M2*xd2-10*sman*xd2-2*tman*xd2+2*uman*xd2+14*M2*xd22-13*sman*xd22+3*tman*xd22-13*uman*xd22-8*M2*xd222+2*sman*xd222+2*tman*xd222-6*uman*xd222-14*M2*xd223+4*sman*xd223+2*tman*xd223-22*M2*xd23+7*sman*xd23+6*tman*xd23-4*M2*xd3+sman*xd3))/(32d0*Pi**2),c0,
     &(GS**4*(2*xc0-8*xc2+2*xc22+2*M2*xd0-2*uman*xd0-2*xd00+8*xd003-2*M2*xd1+2*uman*xd1+17*M2*xd13+sman*xd13-5*uman*xd13+10*M2*xd133-2*sman*xd133+2*uman*xd133-2*M2*xd2+2*uman*xd2+9*M2*xd23+sman*xd23+tman*xd23+2*uman*xd23+8*M2*xd233-2*sman*xd233-2*tman*xd233+6*uman*xd233-5*M2*xd3-3*sman*xd3-tman*xd3+2*uman*xd3-11*M2*xd33+7*sman*xd33+3*tman*xd33+14*M2*xd333-4*sman*xd333-2*tman*xd333))/(32d0*Pi**2),
     &-(GS**4*(2*xc0-8*xc2-2*xc22+7*M2*xd0+sman*xd0+tman*xd0-2*xd00-16*xd001-16*xd002+12*M2*xd1-2*sman*xd1-3*tman*xd1+3*uman*xd1-9*M2*xd11-sman*xd11+2*tman*xd11-uman*xd11-10*M2*xd111+2*sman*xd111-2*uman*xd111-28*M2*xd112+6*sman*xd112+2*tman*xd112-10*uman*xd112-14*M2*xd113+4*sman*xd113+2*tman*xd113-30*M2*xd12-2*sman*xd12+17*tman*xd12-3*uman*xd12-26*M2*xd122+6*sman*xd122+4*tman*xd122-14*uman*xd122-28*M2*xd123+8*sman*xd123+4*tman*xd123-29*M2*xd13+sman*xd13+tman*xd13+6*M2*xd2-2*sman*xd2+2*tman*xd2+4*uman*xd2-21*M2*xd22-sman*xd22+15*tman*xd22-2*uman*xd22-8*M2*xd222+2*sman*xd222+2*tman*xd222-6*uman*xd222-14*M2*xd223+4*sman*xd223+2*tman*xd223-29*M2*xd23+sman*xd23+tman*xd23-7*M2*xd3-sman*xd3-tman*xd3))/(32d0*Pi**2),c0,
     &-(GS**4*(6*xc0-2*xc12-2*xc2+5*M2*xd0+sman*xd0-uman*xd0-10*xd00-4*xd002-4*xd003-5*M2*xd1-sman*xd1+uman*xd1-9*M2*xd12-sman*xd12-3*uman*xd12-10*M2*xd123+2*sman*xd123-2*uman*xd123-5*M2*xd13+2*sman*xd13-3*uman*xd13-2*M2*xd2+2*sman*xd2-2*tman*xd2+4*uman*xd2-M2*xd22-sman*xd22-tman*xd22-10*uman*xd22-8*M2*xd223+2*sman*xd223+2*tman*xd223-6*uman*xd223-2*M2*xd23-5*sman*xd23+7*tman*xd23-5*uman*xd23-14*M2*xd233+4*sman*xd233+2*tman*xd233+6*M2*xd3-3*sman*xd3-3*tman*xd3+uman*xd3-11*M2*xd33+2*sman*xd33+3*tman*xd33))/(32d0*Pi**2),
     &(GS**4*(6*xc0+8*xc1+2*xc12+2*xc2+2*M2*xd0+sman*xd0+2*tman*xd0+6*xd00-4*xd001-12*xd002+3*M2*xd1-3*sman*xd1-2*tman*xd1+3*uman*xd1-5*M2*xd11+2*sman*xd11-3*uman*xd11-10*M2*xd112+2*sman*xd112-2*uman*xd112-5*M2*xd12-3*sman*xd12+4*tman*xd12-3*uman*xd12-18*M2*xd122+4*sman*xd122+2*tman*xd122-8*uman*xd122-14*M2*xd123+4*sman*xd123+2*tman*xd123-11*M2*xd13+2*sman*xd13+3*tman*xd13-2*M2*xd2+8*tman*xd2+6*uman*xd2-10*M2*xd22-5*sman*xd22+11*tman*xd22+3*uman*xd22-8*M2*xd222+2*sman*xd222+2*tman*xd222-6*uman*xd222-14*M2*xd223+4*sman*xd223+2*tman*xd223-14*M2*xd23-sman*xd23-2*tman*xd23-2*M2*xd3-sman*xd3-2*tman*xd3))/(32d0*Pi**2),c0,c0,c0,
     &-(GS**4*(2*xc0+8*xc1-2*xc12-2*xc2+2*M2*xd0-2*uman*xd0+18*xd00-4*xd002-4*xd003-2*M2*xd1+2*uman*xd1-8*M2*xd12+8*uman*xd12-10*M2*xd123+2*sman*xd123-2*uman*xd123-5*M2*xd13-sman*xd13-3*uman*xd13+6*M2*xd2-6*uman*xd2-8*M2*xd22+8*uman*xd22-8*M2*xd223+2*sman*xd223+2*tman*xd223-6*uman*xd223+27*M2*xd23-15*sman*xd23-7*tman*xd23-8*uman*xd23-14*M2*xd233+4*sman*xd233+2*tman*xd233+9*M2*xd3-5*sman*xd3-3*tman*xd3+2*uman*xd3-11*M2*xd33+5*sman*xd33+3*tman*xd33))/(32d0*Pi**2),c0,
     &(GS**4*(-6*xc0+2*xc22+2*M2*xd0-2*uman*xd0-14*xd00-4*xd001-4*xd002-8*xd003+6*M2*xd1-6*uman*xd1-8*M2*xd11+8*uman*xd11-10*M2*xd113+2*sman*xd113-2*uman*xd113-16*M2*xd12+16*uman*xd12-18*M2*xd123+4*sman*xd123+2*tman*xd123-8*uman*xd123+13*M2*xd13-7*sman*xd13-6*tman*xd13+5*uman*xd13-14*M2*xd133+4*sman*xd133+2*tman*xd133+6*M2*xd2-6*uman*xd2-8*M2*xd22+8*uman*xd22-8*M2*xd223+2*sman*xd223+2*tman*xd223-6*uman*xd223+3*M2*xd23-7*sman*xd23+tman*xd23+8*uman*xd23-14*M2*xd233+4*sman*xd233+2*tman*xd233+M2*xd3+3*sman*xd3+5*tman*xd3+2*uman*xd3-3*M2*xd33-3*sman*xd33-5*tman*xd33))/(32d0*Pi**2),
     &(GS**4*(10*xc0+2*xc12+2*xc2+7*M2*xd0+sman*xd0+tman*xd0-22*xd00-4*xd001-12*xd002-2*M2*xd1-tman*xd1+3*uman*xd1-5*M2*xd11-sman*xd11-3*uman*xd11-10*M2*xd112+2*sman*xd112-2*uman*xd112-10*sman*xd12+tman*xd12-17*uman*xd12-18*M2*xd122+4*sman*xd122+2*tman*xd122-8*uman*xd122-14*M2*xd123+4*sman*xd123+2*tman*xd123-11*M2*xd13+5*sman*xd13+3*tman*xd13+14*M2*xd2-10*sman*xd2-6*tman*xd2+4*uman*xd2+3*M2*xd22-9*sman*xd22+7*tman*xd22-18*uman*xd22-8*M2*xd222+2*sman*xd222+2*tman*xd222-6*uman*xd222-14*M2*xd223+4*sman*xd223+2*tman*xd223-37*M2*xd23+9*sman*xd23+9*tman*xd23-7*M2*xd3-sman*xd3-tman*xd3))/(32d0*Pi**2),c0,
     &(GS**4*(14*xc0+2*xc22+11*M2*xd0+sman*xd0+tman*xd0-26*xd00-4*xd001-4*xd002-8*xd003-2*M2*xd1-tman*xd1+3*uman*xd1-9*M2*xd11-sman*xd11-3*uman*xd11-10*M2*xd113+2*sman*xd113-2*uman*xd113-10*M2*xd12-2*sman*xd12-tman*xd12-13*uman*xd12-18*M2*xd123+4*sman*xd123+2*tman*xd123-8*uman*xd123-11*M2*xd13-sman*xd13+5*tman*xd13-6*uman*xd13-14*M2*xd133+4*sman*xd133+2*tman*xd133+6*M2*xd2-4*tman*xd2-2*uman*xd2-M2*xd22-sman*xd22-tman*xd22-10*uman*xd22-8*M2*xd223+2*sman*xd223+2*tman*xd223-6*uman*xd223-13*M2*xd23-sman*xd23+11*tman*xd23-10*uman*xd23-14*M2*xd233+4*sman*xd233+2*tman*xd233+15*M2*xd3-5*sman*xd3-7*tman*xd3-26*M2*xd33+4*sman*xd33+6*tman*xd33))/(32d0*Pi**2)/)      
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo

!     Box 3.1 (21) DOES NOT AGREE WITH FORMCALC EXCEPT FOR Ag (1), (BECAUSE FEYNARTS AMPLITUDE IS WRONG)
      call LTDget(M2,0d0,M2,0d0,tman,uman,M2,0d0,0d0,M2)
      call myCget(0d0,M2,uman,0d0,0d0,M2)

      colors=(/1d0/4d0,0d0,0d0/)
      tempPA=(/(GS**4*(-2*tman*xd00-3*uman*xd00+tman*xd001+uman*xd001+3*tman*xd002+uman*xd002+sman*(xd00+xd002-xd003)-2*uman*xd003+M2*(xd00+2*(xd001+xd003))))/(8d0*Pi**2),
     &(GS**4*(2*xc0+4*xc2+2*M2*xd0-2*tman*xd0-2*xd00-8*xd003-2*M2*xd1+2*tman*xd1-2*M2*xd13+5*tman*xd13+uman*xd13+4*M2*xd133+2*tman*xd133+2*uman*xd133-2*M2*xd2+2*tman*xd2-4*M2*xd23+sman*xd23+7*tman*xd23+uman*xd23+2*sman*xd233+6*tman*xd233+2*uman*xd233+5*M2*xd3-sman*xd3-6*tman*xd3-3*uman*xd3+4*M2*xd33-3*sman*xd33-4*tman*xd33-8*uman*xd33+4*M2*xd333-2*sman*xd333-4*uman*xd333))/(16d0*Pi**2),
     &(GS**4*(-4*xc1-2*xd00+8*xd002+3*M2*xd12-tman*xd12+2*uman*xd12+4*M2*xd122+2*tman*xd122+2*uman*xd122-3*M2*xd2+tman*xd2-2*uman*xd2+3*M2*xd22+4*sman*xd22-3*tman*xd22-4*uman*xd22+2*sman*xd222+6*tman*xd222+2*uman*xd222+4*M2*xd223-2*sman*xd223-4*uman*xd223+2*M2*xd23+sman*xd23-2*uman*xd23))/(16d0*Pi**2),c0,
     &(GS**4*(4*xc0+4*xc2+4*xd00+8*xd001+8*xd002-3*M2*xd1-sman*xd1-uman*xd1+M2*xd11+sman*xd11-tman*xd11+2*M2*xd111+tman*xd111+uman*xd111+4*M2*xd112+sman*xd112+5*tman*xd112+3*uman*xd112+2*M2*xd113-sman*xd113-2*uman*xd113+5*sman*xd12+2*M2*xd122+2*sman*xd122+7*tman*xd122+3*uman*xd122+4*M2*xd123-2*sman*xd123-4*uman*xd123+2*M2*xd13-sman*xd13-2*uman*xd13-3*M2*xd2-sman*xd2-uman*xd2-M2*xd22+4*sman*xd22+tman*xd22+sman*xd222+3*tman*xd222+uman*xd222+2*M2*xd223-sman*xd223-2*uman*xd223+2*M2*xd23-sman*xd23-2*uman*xd23))/(8d0*Pi**2),
     &(GS**4*(4*xc0+8*xc2+3*M2*xd0-tman*xd0+2*uman*xd0+4*xd00+8*xd002-8*xd003-3*M2*xd1+tman*xd1-2*uman*xd1-4*M2*xd12-2*tman*xd12-2*uman*xd12-8*M2*xd123-4*tman*xd123-4*uman*xd123-6*M2*xd13+2*tman*xd13-4*uman*xd13-3*M2*xd2-4*sman*xd2+3*tman*xd2+4*uman*xd2-2*sman*xd22-6*tman*xd22-2*uman*xd22-4*sman*xd223-12*tman*xd223-4*uman*xd223-10*M2*xd23-6*sman*xd23+6*tman*xd23+12*uman*xd23-8*M2*xd233+4*sman*xd233+8*uman*xd233+4*M2*xd3-sman*xd3-2*tman*xd3+6*uman*xd3-4*M2*xd33-2*sman*xd33+4*uman*xd33))/(32d0*Pi**2),c0,
     &-(GS**4*(4*xc0+8*xc2+3*M2*xd0+sman*xd0+uman*xd0+8*xd00+4*xd001+4*xd002-8*xd003-M2*xd1-sman*xd1+tman*xd1-2*M2*xd11-tman*xd11-uman*xd11-4*M2*xd113-2*tman*xd113-2*uman*xd113-2*M2*xd12-sman*xd12-4*tman*xd12-2*uman*xd12-4*M2*xd123-2*sman*xd123-8*tman*xd123-4*uman*xd123-4*M2*xd13-sman*xd13+2*tman*xd13+2*uman*xd13-4*M2*xd133+2*sman*xd133+4*uman*xd133+M2*xd2-4*sman*xd2-tman*xd2-sman*xd22-3*tman*xd22-uman*xd22-2*sman*xd223-6*tman*xd223-2*uman*xd223-7*sman*xd23-2*tman*xd23+2*uman*xd23-4*M2*xd233+2*sman*xd233+4*uman*xd233+4*M2*xd3+3*sman*xd3+4*uman*xd3-4*M2*xd33+2*sman*xd33+4*uman*xd33))/(16d0*Pi**2),
     &-(GS**4*(2*xc1+4*xd00-2*xd002+2*xd003-2*M2*xd12+2*tman*xd12+2*M2*xd123+tman*xd123+uman*xd123+2*M2*xd2-2*tman*xd2-2*M2*xd22+2*tman*xd22+sman*xd223+3*tman*xd223+uman*xd223+M2*xd23-sman*xd23-2*tman*xd23-3*uman*xd23+2*M2*xd233-sman*xd233-2*uman*xd233))/(8d0*Pi**2),c0,
     &-(GS**4*(-4*xc1+2*xd00+2*xd001+6*xd002+2*M2*xd112+tman*xd112+uman*xd112+M2*xd12+sman*xd12-tman*xd12+2*M2*xd122+sman*xd122+4*tman*xd122+2*uman*xd122+2*M2*xd123-sman*xd123-2*uman*xd123-3*M2*xd2-sman*xd2-uman*xd2-M2*xd22+4*sman*xd22+tman*xd22+sman*xd222+3*tman*xd222+uman*xd222+2*M2*xd223-sman*xd223-2*uman*xd223+2*M2*xd23-sman*xd23-2*uman*xd23))/(8d0*Pi**2),c0,c0,c0,
     &(GS**4*(-2*xc0-2*xc2-2*xd001-2*xd002+4*xd003+2*M2*xd1-2*tman*xd1-2*M2*xd11+2*tman*xd11+2*M2*xd113+tman*xd113+uman*xd113-4*M2*xd12+4*tman*xd12+2*M2*xd123+sman*xd123+4*tman*xd123+2*uman*xd123+M2*xd13-sman*xd13-2*tman*xd13-3*uman*xd13+2*M2*xd133-sman*xd133-2*uman*xd133+2*M2*xd2-2*tman*xd2-2*M2*xd22+2*tman*xd22+sman*xd223+3*tman*xd223+uman*xd223+M2*xd23-sman*xd23-2*tman*xd23-3*uman*xd23+2*M2*xd233-sman*xd233-2*uman*xd233))/(8d0*Pi**2),
     &-(GS**4*(4*xc0+4*xc2+4*xd00+4*xd001+12*xd002-3*M2*xd1+tman*xd1-2*uman*xd1+3*M2*xd11-tman*xd11+2*uman*xd11+4*M2*xd112+2*tman*xd112+2*uman*xd112+6*M2*xd12+4*sman*xd12-4*tman*xd12-2*uman*xd12+4*M2*xd122+2*sman*xd122+8*tman*xd122+4*uman*xd122+4*M2*xd123-2*sman*xd123-4*uman*xd123+2*M2*xd13+sman*xd13-2*uman*xd13-3*M2*xd2+tman*xd2-2*uman*xd2+3*M2*xd22+4*sman*xd22-3*tman*xd22-4*uman*xd22+2*sman*xd222+6*tman*xd222+2*uman*xd222+4*M2*xd223-2*sman*xd223-4*uman*xd223+2*M2*xd23+sman*xd23-2*uman*xd23))/(16d0*Pi**2),c0/)      
      
      do n=1,17
      PA_A4(1,n)=PA_A4(1,n)+colors(1)*tempPA(n)
      enddo

!     Box 3.2 (22) AGREES WITH FORMCALC
      call LTDget(M2,0d0,M2,0d0,tman,uman,0d0,M2,M2,0d0)
      call myCget(0d0,M2,uman,M2,M2,0d0)


      colors=(/1d0/4d0,0d0,0d0/)
      tempPA=(/-(GS**4*(4*tman*xd00+3*uman*xd00+tman*xd001+uman*xd001+3*tman*xd002+uman*xd002-2*uman*xd003-sman*(xd00-xd002+xd003)+M2*(xd00+2*(xd001+xd003))))/(8d0*Pi**2),
     &-(GS**4*(4*xc2+2*xd00-8*xd003+11*M2*xd13+2*tman*xd13-uman*xd13+4*M2*xd133+2*tman*xd133+2*uman*xd133+5*M2*xd23-sman*xd23+8*tman*xd23-uman*xd23+2*sman*xd233+6*tman*xd233+2*uman*xd233+6*M2*xd3+4*tman*xd3-2*uman*xd3+8*M2*xd33-4*sman*xd33+8*tman*xd33+4*M2*xd333-2*sman*xd333-4*uman*xd333))/(16d0*Pi**2),
     &(GS**4*(2*xc0+4*xc1+4*M2*xd0-4*uman*xd0-2*xd00-8*xd002+6*M2*xd1-2*uman*xd1+10*M2*xd12-tman*xd12-5*uman*xd12-4*M2*xd122-2*tman*xd122-2*uman*xd122+11*M2*xd2-sman*xd2-2*tman*xd2-13*uman*xd2+6*M2*xd22-3*sman*xd22-7*tman*xd22-11*uman*xd22-2*sman*xd222-6*tman*xd222-2*uman*xd222-4*M2*xd223+2*sman*xd223+4*uman*xd223+2*M2*xd23+sman*xd23-2*uman*xd23+2*M2*xd3-2*uman*xd3))/(16d0*Pi**2),c0,
     &-(GS**4*(4*xc2-2*sman*xd0+2*uman*xd0-4*xd00+8*xd001+8*xd002-3*sman*xd1+2*tman*xd1+5*uman*xd1-sman*xd11+3*tman*xd11+4*uman*xd11+tman*xd111+uman*xd111+sman*xd112+5*tman*xd112+3*uman*xd112-sman*xd113-2*uman*xd113-sman*xd12+4*tman*xd12+8*uman*xd12+2*sman*xd122+7*tman*xd122+3*uman*xd122-2*sman*xd123-4*uman*xd123+sman*xd13+2*uman*xd13-3*sman*xd2-2*tman*xd2+5*uman*xd2+tman*xd22+4*uman*xd22+sman*xd222+3*tman*xd222+uman*xd222-sman*xd223-2*uman*xd223+sman*xd23+2*uman*xd23+2*sman*xd3+4*uman*xd3-M2*(2*xd0+7*xd1+3*xd11-2*xd111-4*xd112-2*xd113+4*xd12-2*xd122-4*xd123+2*xd13+3*xd2+xd22-2*xd223+2*xd23+4*xd3)))/(8d0*Pi**2),
     &(GS**4*(4*xc0+8*xc1+6*M2*xd0+4*tman*xd0-2*uman*xd0+4*xd00-8*xd002+8*xd003+11*M2*xd1+2*tman*xd1-uman*xd1+22*M2*xd12+4*tman*xd12-2*uman*xd12+8*M2*xd123+4*tman*xd123+4*uman*xd123+4*M2*xd13+2*tman*xd13+2*uman*xd13+17*M2*xd2-sman*xd2+16*tman*xd2-5*uman*xd2+10*M2*xd22-2*sman*xd22+16*tman*xd22-2*uman*xd22+4*sman*xd223+12*tman*xd223+4*uman*xd223+16*M2*xd23-6*sman*xd23+22*tman*xd23+2*uman*xd23+8*M2*xd233-4*sman*xd233-8*uman*xd233+8*M2*xd3-4*sman*xd3+8*tman*xd3+4*M2*xd33-2*sman*xd33-4*uman*xd33))/(32d0*Pi**2),c0,
     &(GS**4*(4*xc2-4*tman*xd0+2*uman*xd0-8*xd00+4*xd001+4*xd002-8*xd003-6*tman*xd1+3*uman*xd1-2*tman*xd11+uman*xd11-2*tman*xd113-2*uman*xd113+sman*xd12-10*tman*xd12+2*uman*xd12-2*sman*xd123-8*tman*xd123-4*uman*xd123+4*sman*xd13-10*tman*xd13-2*uman*xd13+2*sman*xd133+4*uman*xd133+sman*xd2-12*tman*xd2+3*uman*xd2+sman*xd22-8*tman*xd22+uman*xd22-2*sman*xd223-6*tman*xd223-2*uman*xd223+2*sman*xd23-14*tman*xd23-2*uman*xd23+2*sman*xd233+4*uman*xd233+4*sman*xd3-8*tman*xd3+2*sman*xd33+4*uman*xd33-M2*(6*xd0+17*xd1+11*xd11+4*xd113+16*xd12+4*xd123+12*xd13+4*xd133+11*xd2+5*xd22+8*xd23+4*xd233+8*xd3+4*xd33)))/(16d0*Pi**2),
     &(GS**4*(-2*xc2-4*xd00-2*xd002+2*xd003+2*M2*xd123+tman*xd123+uman*xd123-6*M2*xd13+2*uman*xd13+sman*xd223+3*tman*xd223+uman*xd223-3*M2*xd23+sman*xd23+2*tman*xd23+5*uman*xd23+2*M2*xd233-sman*xd233-2*uman*xd233-4*M2*xd3+4*uman*xd3-2*M2*xd33+2*uman*xd33))/(8d0*Pi**2),c0,
     &(GS**4*(2*xc2+4*uman*xd0-2*xd00+2*xd001+6*xd002+6*uman*xd1+2*uman*xd11+tman*xd112+uman*xd112+sman*xd12+3*tman*xd12+8*uman*xd12+sman*xd122+4*tman*xd122+2*uman*xd122-sman*xd123-2*uman*xd123+2*uman*xd13+sman*xd2+2*tman*xd2+9*uman*xd2+2*sman*xd22+5*tman*xd22+6*uman*xd22+sman*xd222+3*tman*xd222+uman*xd222-sman*xd223-2*uman*xd223-sman*xd23+2*uman*xd3-M2*(4*xd0+10*xd1+6*xd11-2*xd112+7*xd12-2*xd122-2*xd123+2*xd13+7*xd2+3*xd22-2*xd223+2*xd3)))/(8d0*Pi**2),c0,c0,c0,
     &-(GS**4*(-4*xc2-2*xd001-2*xd002+4*xd003+2*M2*xd113+tman*xd113+uman*xd113+2*M2*xd123+sman*xd123+4*tman*xd123+2*uman*xd123-5*M2*xd13-sman*xd13+2*tman*xd13+3*uman*xd13+2*M2*xd133-sman*xd133-2*uman*xd133+sman*xd223+3*tman*xd223+uman*xd223-M2*xd23-sman*xd23-2*tman*xd23+3*uman*xd23+2*M2*xd233-sman*xd233-2*uman*xd233-2*M2*xd3-2*sman*xd3+2*uman*xd3-4*M2*xd33+2*sman*xd33+4*uman*xd33))/(8d0*Pi**2),
     &-(GS**4*(4*xc0+8*xc1+2*M2*xd0+2*sman*xd0-2*uman*xd0+4*xd00-4*xd001-12*xd002+5*M2*xd1+sman*xd1-2*tman*xd1-3*uman*xd1-2*M2*xd11-tman*xd11-uman*xd11-4*M2*xd112-2*tman*xd112-2*uman*xd112+8*M2*xd12+sman*xd12-8*tman*xd12-8*uman*xd12-4*M2*xd122-2*sman*xd122-8*tman*xd122-4*uman*xd122-4*M2*xd123+2*sman*xd123+4*uman*xd123-2*M2*xd13+sman*xd13+2*uman*xd13+5*M2*xd2+5*sman*xd2+2*tman*xd2-7*uman*xd2+2*M2*xd22+sman*xd22+tman*xd22-7*uman*xd22-2*sman*xd222-6*tman*xd222-2*uman*xd222-4*M2*xd223+2*sman*xd223+4*uman*xd223+6*M2*xd23-3*sman*xd23-6*uman*xd23+4*M2*xd3-2*sman*xd3-4*uman*xd3))/(16d0*Pi**2),c0/)      
      
      do n=1,17
      PA_A4(1,n)=PA_A4(1,n)+colors(1)*tempPA(n)
      enddo

      endif

!     SquarkExchangeRight (16) 
      do k=1,2 
!     First part: Sum over all possible squarks running in the loop
      do gammaType=3,4
      do gammaGen=1,3

      call myCget(0d0,sman,0d0,MSf2(k,gammaType,gammaGen),MSf2(k,gammaType,gammaGen),MSf2(k,gammaType,gammaGen))

      tempPA=(/(GS**4)/(2d0*Pi**2)*xc00, (GS**4)/(4d0*Pi**2)*(xc1+2*xc11),(GS**4)/(4*Pi**2)*(xc1+2*xc11),c0,c0,(-GS**4)/(8d0*Pi**2)*(xc0+4*(xc1+xc12)),
     &c0,c0,(-GS**4)/(2d0*Pi**2)*xc12,c0,c0,c0,c0,c0,c0,c0,c0/)


      tempMixing=(USf(j,1,type1,gen1)*USf(i,1,type1,gen1)-USf(j,2,type1,gen1)*USf(i,2,type1,gen1))*
     &   (USf(k,1,gammaType,gammaGen)* USf(k,1,gammaType,gammaGen)- USf(k,2,gammaType,gammaGen)*USf(k,2,gammaType,gammaGen))

      do n=1,9
      PA_A4(2,n)=PA_A4(2,n)+tempPA(n)*tempMixing/4d0
      enddo

      enddo
      enddo 

!     Second part: sum over the two mass eigenstates which belong to the same generation and type of the incoming squarks
      call myCget(0d0,sman,0d0,MSf2(k,type1,gen1), MSf2(k,type1,gen1), MSf2(k,type1,gen1))

      tempPA=(/(GS**4)/(2d0*Pi**2)*xc00, (GS**4)/(4d0*Pi**2)*(xc1+2*xc11),(GS**4)/(4d0*Pi**2)*(xc1+2*xc11),c0,c0,(-GS**4)/(8d0*Pi**2)*(xc0+4*(xc1+xc12)),
     &c0,c0,(-GS**4)/(2d0*Pi**2)*xc12,c0,c0,c0,c0,c0,c0,c0,c0/)

      tempMixing=(USf(j,1,type1,gen1)*USf(k,1,type1,gen1)-USf(j,2,type1,gen1)*USf(k,2,type1,gen1))*
     &   (USf(k,1,type1,gen1)* USf(i,1,type1,gen1)- USf(k,2,type1,gen1)*USf(i,2,type1,gen1))

      colors=(/CF/(2d0*Nc),-1d0/(4d0*Nc),0d0/)

      do m=1,2
      do n=1,9
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)*tempMixing
      enddo
      enddo
      
      enddo


!     SquarkLoop (1) 
      do k=1,2 
!     First part: Sum over all possible squarks running in the loop
      do gammaType=3,4
      do gammaGen=1,3

      call myBget(sman,MSf2(k,gammaType,gammaGen),MSf2(k,gammaType,gammaGen))

      temp=-GS**4/(16d0*Pi**2)*xb0

      tempMixing=(USf(j,1,type1,gen1)*USf(i,1,type1,gen1)-USf(j,2,type1,gen1)*USf(i,2,type1,gen1))*
     &   (USf(k,1,gammaType,gammaGen)* USf(k,1,gammaType,gammaGen)- USf(k,2,gammaType,gammaGen)*USf(k,2,gammaType,gammaGen))

      PA_A4(2,1)=PA_A4(2,1)+temp*tempMixing/2d0

      enddo
      enddo

!     Second part: sum over the two mass eigenstates which belong to the same generation and type of the incoming squarks
      call myBget(sman,MSf2(k,type1,gen1),MSf2(k,type1,gen1))

      temp=-GS**4/(16d0*Pi**2)*xb0

      tempMixing=(USf(j,1,type1,gen1)*USf(k,1,type1,gen1)-USf(j,2,type1,gen1)*USf(k,2,type1,gen1))*
     &   (USf(k,1,type1,gen1)* USf(i,1,type1,gen1)- USf(k,2,type1,gen1)*USf(i,2,type1,gen1))

      colors=(/CF/Nc,-1d0/(2d0*Nc),0d0/)

      do m=1,2
      PA_A4(m,1)=PA_A4(m,1)+colors(m)*temp*tempMixing
      enddo     
      enddo


!     Box 4 (23) A43 (17) DOES NOT AGREE WITH FORMCALC
      call LTDget(M2j,0d0,0d0,M2,tman,sman,MQ2,MGl2,MGl2,MGl2)
      call myCget(0d0,0d0,sman,MGl2,MGl2,MGl2)
      call myBget(0d0,MGl2,MGl2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/1d0/2d0,Nc/4d0,-Nc/4d0/)
      tempPA=(/-(GS**4*(MGl*MQ*tempMixing*(2*xc0+M2*xd0+M2j*xd0-2*MGl2*xd0+2*MQ2*xd0+tman*xd0-uman*xd0-8*xd00+2*M2j*xd1+2*tman*xd1+4*tman*xd2+4*M2*xd3+2*M2j*xd3-2*sman*xd3-2*uman*xd3)-(2*xb0+2*MQ2*xc0+sman*xc0-uman*xc0-4*xc00+M2*xc1-uman*xc1+2*M2*xc2+sman*xc2-tman*xc2-uman*xc2+2*MQ**4*xd0+M2*MQ2*xd0-2*MGl2*MQ2*xd0+MQ2*tman*xd0-MQ2*uman*xd0
     &+4*M2*xd00-4*MGl2*xd00-4*MQ2*xd00-4*sman*xd00-M2*MGl2*xd1+M2*MQ2*xd1+MGl2*sman*xd1-MQ2*sman*xd1+MGl2*tman*xd1+MQ2*tman*xd1+M2*tman*xd12-sman*tman*xd12+M2**2*xd13-2*M2*sman*xd13+sman**2*xd13-2*M2*MGl2*xd2+2*M2*MQ2*xd2+MGl2*sman*xd2-MQ2*sman*xd2+M2*tman*xd2+MGl2*tman*xd2+3*MQ2*tman*xd2+MGl2*uman*xd2-MQ2*uman*xd2-tman*uman*xd2+2*M2*tman*xd22
     &-sman*tman*xd22-tman*uman*xd22+2*M2**2*xd23-3*M2*sman*xd23+sman**2*xd23+2*M2*tman*xd23-M2*uman*xd23+sman*uman*xd23+M2j**2*(xd1+2*xd11+3*xd12+xd13+xd2+xd22+xd23)-M2*MGl2*xd3+5*M2*MQ2*xd3-2*MQ2*sman*xd3+M2*tman*xd3-MGl2*uman*xd3-MQ2*uman*xd3+2*M2**2*xd33-2*M2*sman*xd33-M2j*(xc0-4*xd00-2*M2*xd1+2*MGl2*xd1+sman*xd1+uman*xd1-2*M2*xd11+2*sman*xd11
     &-5*M2*xd12+3*sman*xd12-tman*xd12+2*uman*xd12-6*M2*xd13+2*sman*xd13-2*M2*xd2+2*MGl2*xd2+sman*xd2+uman*xd2-2*M2*xd22+sman*xd22-tman*xd22+uman*xd22-5*M2*xd23+2*sman*xd23+uman*xd23-M2*xd3-MQ2*(xd0+2*(2*xd1+xd2+xd3))-2*M2*xd33))*Krondelta(i,j)))/(8d0*Pi**2),c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd2+xd22)+(2*xc1+2*xc11+(-M2-M2j+2*MGl2+2*MQ2+sman)*(xd2+xd22))*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)+(2*xc2+2*xc22+2*M2j*xd13+M2j*xd23+tman*xd23+2*MQ2*xd3+2*MGl2*xd33+2*MQ2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+3*xd1+2*xd11+4*xd12+3*xd2+2*xd22)+(2*xc22+2*MGl2*xd1+2*MGl2*xd11+M2*xd12+4*MGl2*xd12-uman*xd12+2*M2*xd13+M2*xd2+2*MGl2*xd2-uman*xd2+M2*xd22+2*MGl2*xd22-uman*xd22+2*MQ2*(xd0+2*xd1+xd11+2*xd12+2*xd2+xd22)+2*M2*xd23+2*M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+4*xd23+2*xd3)-(xc0+xc1+2*xc12+xc2+MQ2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2-MGl2*xd2+MQ2*xd2+tman*xd2+M2j*xd22+tman*xd22+2*MGl2*xd23+2*MQ2*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*xd12+6*xd2+4*xd22))+(-2*xc12-xc2+MQ2*xd0+MGl2*xd1+MQ2*xd1-M2*xd12-M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+sman*xd12-M2*xd2-M2j*xd2+3*MGl2*xd2+3*MQ2*xd2+sman*xd2-M2*xd22-M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3))-(xc0+xc1+2*xc12+xc2+MQ2*xd0+M2j*xd1+MGl2*xd2+MQ2*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0-xc2-2*xc22+2*MQ2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1-M2*xd13-M2j*xd13+2*MGl2*xd13+2*MQ2*xd13+sman*xd13+MGl2*xd2+MQ2*xd2+tman*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+2*MGl2*xd3+2*MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*(xd12+xd2+xd22)))+(-2*xc12-xc2+MQ2*xd0+MGl2*xd1+MQ2*xd1+2*MGl2*xd12+2*MQ2*xd12+M2*xd2+MGl2*xd2+3*MQ2*xd2-uman*xd2+M2*xd22+2*MGl2*xd22+2*MQ2*xd22-uman*xd22+2*M2*xd23+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(2*xd13+2*xd23+xd3)+(-xc0-xc2-2*xc22+4*xd00+M2j*xd1-MGl2*xd1+MQ2*xd1+2*M2j*xd11+3*M2j*xd12+tman*xd12+M2*xd13+M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-sman*xd13+M2j*xd2-MGl2*xd2+MQ2*xd2+M2j*xd22+tman*xd22+2*M2*xd23+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-sman*xd23-uman*xd23+2*MQ2*xd3+2*M2*xd33)*Krondelta(i,j)))/(4d0*Pi**2)/)     
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo


!     Box 4 u (24) A43 (17) DOES NOT AGREE WITH FORMCALC
      call LTDget(M2j,0d0,0d0,M2,uman,sman,MQ2,MGl2,MGl2,MGl2)
      call myCget(0d0,0d0,sman,MGl2,MGl2,MGl2)
      call myBget(0d0,MGl2,MGl2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/1d0/2d0,Nc/4d0,Nc/4d0/)
      tempPA=(/-(GS**4*(MGl*MQ*tempMixing*(2*xc0+M2*xd0+M2j*xd0-2*MGl2*xd0+2*MQ2*xd0-tman*xd0+uman*xd0-8*xd00+2*M2j*xd1+2*uman*xd1+4*uman*xd2+4*M2*xd3+2*M2j*xd3-2*sman*xd3-2*tman*xd3)-(2*xb0+2*MQ2*xc0+sman*xc0-tman*xc0-4*xc00+M2*xc1-tman*xc1+2*M2*xc2+sman*xc2-tman*xc2-uman*xc2+2*MQ**4*xd0+M2*MQ2*xd0-2*MGl2*MQ2*xd0-MQ2*tman*xd0+MQ2*uman*xd0
     &+4*M2*xd00-4*MGl2*xd00-4*MQ2*xd00-4*sman*xd00-M2*MGl2*xd1+M2*MQ2*xd1+MGl2*sman*xd1-MQ2*sman*xd1+MGl2*uman*xd1+MQ2*uman*xd1+M2*uman*xd12-sman*uman*xd12+M2**2*xd13-2*M2*sman*xd13+sman**2*xd13-2*M2*MGl2*xd2+2*M2*MQ2*xd2+MGl2*sman*xd2-MQ2*sman*xd2+MGl2*tman*xd2-MQ2*tman*xd2+M2*uman*xd2+MGl2*uman*xd2+3*MQ2*uman*xd2-tman*uman*xd2+2*M2*uman*xd22
     &-sman*uman*xd22-tman*uman*xd22+2*M2**2*xd23-3*M2*sman*xd23+sman**2*xd23-M2*tman*xd23+sman*tman*xd23+2*M2*uman*xd23+M2j**2*(xd1+2*xd11+3*xd12+xd13+xd2+xd22+xd23)-M2*MGl2*xd3+5*M2*MQ2*xd3-2*MQ2*sman*xd3-MGl2*tman*xd3-MQ2*tman*xd3+M2*uman*xd3+2*M2**2*xd33-2*M2*sman*xd33-M2j*(xc0-4*xd00-2*M2*xd1+2*MGl2*xd1+sman*xd1+tman*xd1-2*M2*xd11+2*sman*xd11
     &-5*M2*xd12+3*sman*xd12+2*tman*xd12-uman*xd12-6*M2*xd13+2*sman*xd13-2*M2*xd2+2*MGl2*xd2+sman*xd2+tman*xd2-2*M2*xd22+sman*xd22+tman*xd22-uman*xd22-5*M2*xd23+2*sman*xd23+tman*xd23-M2*xd3-MQ2*(xd0+2*(2*xd1+xd2+xd3))-2*M2*xd33))*Krondelta(i,j)))/(8d0*Pi**2),
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd2+xd22)+(2*xc1+2*xc11+(-M2-M2j+2*MGl2+2*MQ2+sman)*(xd2+xd22))*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)+(2*xc2+2*xc22+2*M2j*xd13+M2j*xd23+uman*xd23+2*MQ2*xd3+2*MGl2*xd33+2*MQ2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+3*xd1+2*xd11+4*xd12+3*xd2+2*xd22)+(2*xc22+2*MGl2*xd1+2*MGl2*xd11+M2*xd12+4*MGl2*xd12-tman*xd12+2*M2*xd13+M2*xd2+2*MGl2*xd2-tman*xd2+M2*xd22+2*MGl2*xd22-tman*xd22+2*MQ2*(xd0+2*xd1+xd11+2*xd12+2*xd2+xd22)+2*M2*xd23+2*M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3))-(xc0+xc1+2*xc12+xc2+MQ2*xd0+M2j*xd1+MGl2*xd2+MQ2*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*(xd12+xd2+xd22)))+(-2*xc12-xc2+MQ2*xd0+MGl2*xd1+MQ2*xd1+2*MGl2*xd12+2*MQ2*xd12+M2*xd2+MGl2*xd2+3*MQ2*xd2-tman*xd2+M2*xd22+2*MGl2*xd22+2*MQ2*xd22-tman*xd22+2*M2*xd23+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+4*xd23+2*xd3)-(xc0+xc1+2*xc12+xc2+MQ2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2-MGl2*xd2+MQ2*xd2+uman*xd2+M2j*xd22+uman*xd22+2*MGl2*xd23+2*MQ2*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-2*MGl*MQ*tempMixing*(2*xd13+2*xd23+xd3)+(-xc0-xc2-2*xc22+4*xd00+M2j*xd1-MGl2*xd1+MQ2*xd1+2*M2j*xd11+3*M2j*xd12+uman*xd12+M2*xd13+M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-sman*xd13+M2j*xd2-MGl2*xd2+MQ2*xd2+M2j*xd22+uman*xd22+2*M2*xd23+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-sman*xd23-tman*xd23+2*MQ2*xd3+2*M2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*xd12+6*xd2+4*xd22))+(-2*xc12-xc2+MQ2*xd0+MGl2*xd1+MQ2*xd1-M2*xd12-M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+sman*xd12-M2*xd2-M2j*xd2+3*MGl2*xd2+3*MQ2*xd2+sman*xd2-M2*xd22-M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0-xc2-2*xc22+2*MQ2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1-M2*xd13-M2j*xd13+2*MGl2*xd13+2*MQ2*xd13+sman*xd13+MGl2*xd2+MQ2*xd2+uman*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+2*MGl2*xd3+2*MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2)/)     
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo
  

!     Box 5 (25) A33 (4), A23 (10) AND A43 (17) DO NOT AGREE WITH FORMCALC
      call LTDget(M2j,0d0,0d0,M2,tman,sman,MGl2,MQ2,MQ2,MQ2)
      call myCget(0d0,0d0,sman,MQ2,MQ2,MQ2)
      call myBget(0d0,MQ2,MQ2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/CF/(2d0*Nc),-1d0/(4d0*Nc),1d0/(4d0*Nc)/)
      tempPA=(/(GS**4*(-(MGl*MQ*tempMixing*(2*xc0+M2*xd0+M2j*xd0+2*MGl2*xd0-2*MQ2*xd0+tman*xd0-uman*xd0-8*xd00+2*M2j*xd1+2*tman*xd1+4*tman*xd2+4*M2*xd3+2*M2j*xd3-2*sman*xd3-2*uman*xd3))+(2*xb0-M2j*xc0+sman*xc0-uman*xc0-4*xc00+M2*xc1-uman*xc1+2*M2*xc2+sman*xc2-tman*xc2-uman*xc2+2*MGl2**2*xd0+4*M2*xd00+4*M2j*xd00-4*MQ2*xd00-4*sman*xd00+2*M2*M2j*xd1
     &+M2j**2*xd1-M2*MQ2*xd1-2*M2j*MQ2*xd1-M2j*sman*xd1+MQ2*sman*xd1+MQ2*tman*xd1-M2j*uman*xd1+2*M2*M2j*xd11+2*M2j**2*xd11-2*M2j*sman*xd11+5*M2*M2j*xd12+3*M2j**2*xd12-3*M2j*sman*xd12+M2*tman*xd12+M2j*tman*xd12-sman*tman*xd12-2*M2j*uman*xd12+M2**2*xd13+6*M2*M2j*xd13+M2j**2*xd13-2*M2*sman*xd13-2*M2j*sman*xd13+sman**2*xd13+2*M2*M2j*xd2+M2j**2*xd2-2*M2*MQ2*xd2
     &-2*M2j*MQ2*xd2-M2j*sman*xd2+MQ2*sman*xd2+M2*tman*xd2+MQ2*tman*xd2-M2j*uman*xd2+MQ2*uman*xd2-tman*uman*xd2+2*M2*M2j*xd22+M2j**2*xd22-M2j*sman*xd22+2*M2*tman*xd22+M2j*tman*xd22-sman*tman*xd22-M2j*uman*xd22-tman*uman*xd22+2*M2**2*xd23+5*M2*M2j*xd23+M2j**2*xd23-3*M2*sman*xd23-2*M2j*sman*xd23+sman**2*xd23+2*M2*tman*xd23-M2*uman*xd23-M2j*uman*xd23
     &+sman*uman*xd23+M2*M2j*xd3-M2*MQ2*xd3+M2*tman*xd3-MQ2*uman*xd3+MGl2*(2*xc0+M2j*xd0-2*MQ2*xd0+tman*xd0-uman*xd0-4*xd00+4*M2j*xd1-sman*xd1+tman*xd1+2*M2j*xd2-sman*xd2+3*tman*xd2-uman*xd2+2*M2j*xd3-2*sman*xd3-uman*xd3+M2*(xd0+xd1+2*xd2+5*xd3))+2*M2**2*xd33+2*M2*M2j*xd33-2*M2*sman*xd33)*Krondelta(i,j)))/(8d0*Pi**2),c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd2+xd22)+(2*xc1+2*xc11+(-M2-M2j+2*MGl2+2*MQ2+sman)*(xd2+xd22))*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)+(2*xc2+2*xc22+2*M2j*xd13+M2j*xd23+tman*xd23+2*MGl2*xd3+2*MGl2*xd33+2*MQ2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+3*xd1+2*xd11+4*xd12+3*xd2+2*xd22)+(2*xc22+2*MQ2*xd1+2*MQ2*xd11+M2*xd12+4*MQ2*xd12-uman*xd12+2*M2*xd13+M2*xd2+2*MQ2*xd2-uman*xd2+M2*xd22+2*MQ2*xd22-uman*xd22+2*MGl2*(xd0+2*xd1+xd11+2*xd12+2*xd2+xd22)+2*M2*xd23+2*M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+4*xd23+2*xd3)-(xc0+xc1+2*xc12+xc2+MGl2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2+MGl2*xd2-MQ2*xd2+tman*xd2+M2j*xd22+tman*xd22+2*MGl2*xd23+2*MQ2*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*xd12+6*xd2+4*xd22))+(-2*xc12-xc2+MGl2*xd0+MGl2*xd1+MQ2*xd1-M2*xd12-M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+sman*xd12-M2*xd2-M2j*xd2+3*MGl2*xd2+3*MQ2*xd2+sman*xd2-M2*xd22-M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3))-(xc0+xc1+2*xc12+xc2+MGl2*xd0+M2j*xd1+MGl2*xd2+MQ2*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0-xc2-2*xc22+2*MGl2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1-M2*xd13-M2j*xd13+2*MGl2*xd13+2*MQ2*xd13+sman*xd13+MGl2*xd2+MQ2*xd2+tman*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+2*MGl2*xd3+2*MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*(xd12+xd2+xd22)))+(-2*xc12-xc2+MGl2*xd0+MGl2*xd1+MQ2*xd1+2*MGl2*xd12+2*MQ2*xd12+M2*xd2+3*MGl2*xd2+MQ2*xd2-uman*xd2+M2*xd22+2*MGl2*xd22+2*MQ2*xd22-uman*xd22+2*M2*xd23+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(2*xd13+2*xd23+xd3)+(-xc0-xc2-2*xc22+4*xd00+M2j*xd1+MGl2*xd1-MQ2*xd1+2*M2j*xd11+3*M2j*xd12+tman*xd12+M2*xd13+M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-sman*xd13+M2j*xd2+MGl2*xd2-MQ2*xd2+M2j*xd22+tman*xd22+2*M2*xd23+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-sman*xd23-uman*xd23+2*MGl2*xd3+2*M2*xd33)*Krondelta(i,j)))/(4d0*Pi**2)/)     
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo


!     Box 5 u (26) A14 (8) AND A34 (14) DO NOT AGREE WITH FORMCALC
      call LTDget(M2j,0d0,0d0,M2,uman,sman,MGl2,MQ2,MQ2,MQ2)
      call myCget(0d0,0d0,sman,MQ2,MQ2,MQ2)
      call myBget(0d0,MQ2,MQ2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/CF/(2d0*Nc),-1d0/(4d0*Nc),-1d0/(4d0*Nc)/)
      tempPA=(/(GS**4*(-(MGl*MQ*tempMixing*(2*xc0+M2*xd0+M2j*xd0+2*MGl2*xd0-2*MQ2*xd0-tman*xd0+uman*xd0-8*xd00+2*M2j*xd1+2*uman*xd1+4*uman*xd2+4*M2*xd3+2*M2j*xd3-2*sman*xd3-2*tman*xd3))+(2*xb0-M2j*xc0+sman*xc0-tman*xc0-4*xc00+M2*xc1-tman*xc1+2*M2*xc2+sman*xc2-tman*xc2-uman*xc2+2*MGl2**2*xd0+4*M2*xd00+4*M2j*xd00-4*MQ2*xd00-4*sman*xd00+2*M2*M2j*xd1
     &+M2j**2*xd1-M2*MQ2*xd1-2*M2j*MQ2*xd1-M2j*sman*xd1+MQ2*sman*xd1-M2j*tman*xd1+MQ2*uman*xd1+2*M2*M2j*xd11+2*M2j**2*xd11-2*M2j*sman*xd11+5*M2*M2j*xd12+3*M2j**2*xd12-3*M2j*sman*xd12-2*M2j*tman*xd12+M2*uman*xd12+M2j*uman*xd12-sman*uman*xd12+M2**2*xd13+6*M2*M2j*xd13+M2j**2*xd13-2*M2*sman*xd13-2*M2j*sman*xd13+sman**2*xd13+2*M2*M2j*xd2+M2j**2*xd2
     &-2*M2*MQ2*xd2-2*M2j*MQ2*xd2-M2j*sman*xd2+MQ2*sman*xd2-M2j*tman*xd2+MQ2*tman*xd2+M2*uman*xd2+MQ2*uman*xd2-tman*uman*xd2+2*M2*M2j*xd22+M2j**2*xd22-M2j*sman*xd22-M2j*tman*xd22+2*M2*uman*xd22+M2j*uman*xd22-sman*uman*xd22-tman*uman*xd22+2*M2**2*xd23+5*M2*M2j*xd23+M2j**2*xd23-3*M2*sman*xd23-2*M2j*sman*xd23+sman**2*xd23-M2*tman*xd23-M2j*tman*xd23
     &+sman*tman*xd23+2*M2*uman*xd23+M2*M2j*xd3-M2*MQ2*xd3-MQ2*tman*xd3+M2*uman*xd3+MGl2*(2*xc0+M2j*xd0-2*MQ2*xd0-tman*xd0+uman*xd0-4*xd00+4*M2j*xd1-sman*xd1+uman*xd1+2*M2j*xd2-sman*xd2-tman*xd2+3*uman*xd2+2*M2j*xd3-2*sman*xd3-tman*xd3+M2*(xd0+xd1+2*xd2+5*xd3))+2*M2**2*xd33+2*M2*M2j*xd33-2*M2*sman*xd33)*Krondelta(i,j)))/(8d0*Pi**2),
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd2+xd22)+(2*xc1+2*xc11+(-M2-M2j+2*MGl2+2*MQ2+sman)*(xd2+xd22))*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)+(2*xc2+2*xc22+2*M2j*xd13+M2j*xd23+uman*xd23+2*MGl2*xd3+2*MGl2*xd33+2*MQ2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+3*xd1+2*xd11+4*xd12+3*xd2+2*xd22)+(2*xc22+2*MQ2*xd1+2*MQ2*xd11+M2*xd12+4*MQ2*xd12-tman*xd12+2*M2*xd13+M2*xd2+2*MQ2*xd2-tman*xd2+M2*xd22+2*MQ2*xd22-tman*xd22+2*MGl2*(xd0+2*xd1+xd11+2*xd12+2*xd2+xd22)+2*M2*xd23+2*M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3))-(xc0+xc1+2*xc12+xc2+MGl2*xd0+M2j*xd1+MGl2*xd2+MQ2*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*(xd12+xd2+xd22)))+(-2*xc12-xc2+MGl2*xd0+MGl2*xd1+MQ2*xd1+2*MGl2*xd12+2*MQ2*xd12+M2*xd2+3*MGl2*xd2+MQ2*xd2-tman*xd2+M2*xd22+2*MGl2*xd22+2*MQ2*xd22-tman*xd22+2*M2*xd23+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,
     &(GS**4*(MGl*MQ*tempMixing*(xd0+4*xd23+2*xd3)-(xc0+xc1+2*xc12+xc2+MGl2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2+MGl2*xd2-MQ2*xd2+uman*xd2+M2j*xd22+uman*xd22+2*MGl2*xd23+2*MQ2*xd23+MGl2*xd3+MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-2*MGl*MQ*tempMixing*(2*xd13+2*xd23+xd3)+(-xc0-xc2-2*xc22+4*xd00+M2j*xd1+MGl2*xd1-MQ2*xd1+2*M2j*xd11+3*M2j*xd12+uman*xd12+M2*xd13+M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-sman*xd13+M2j*xd2+MGl2*xd2-MQ2*xd2+M2j*xd22+uman*xd22+2*M2*xd23+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-sman*xd23-tman*xd23+2*MGl2*xd3+2*M2*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*xd1+4*xd12+6*xd2+4*xd22))+(-2*xc12-xc2+MGl2*xd0+MGl2*xd1+MQ2*xd1-M2*xd12-M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+sman*xd12-M2*xd2-M2j*xd2+3*MGl2*xd2+3*MQ2*xd2+sman*xd2-M2*xd22-M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2*xd3)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0-xc2-2*xc22+2*MGl2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1-M2*xd13-M2j*xd13+2*MGl2*xd13+2*MQ2*xd13+sman*xd13+MGl2*xd2+MQ2*xd2+uman*xd2-M2*xd23-M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23+2*MGl2*xd3+2*MQ2*xd3)*Krondelta(i,j)))/(4d0*Pi**2)/)      
      
      do n=1,17
      do m=1,3
      PA_A4(m,n)=PA_A4(m,n)+colors(m)*tempPA(n)
      enddo
      enddo


!     Box 6.1 (27) DOES NOT AGREE WITH FORMCALC EXCEPT FOR Ag (1), (BECAUSE FEYNARTS AMPLITUDE IS WRONG)
      call LTDget(M2j,0d0,M2,0d0,tman,uman,MQ2,MGl2,MGl2,MQ2)
      call myCget(0d0,M2,uman,MGl2,MGl2,MQ2)     
      call myBget(M2,MGl2,MQ2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/1d0/4d0,0d0,0d0/)
      tempPA=(/-(GS**4*(MGl*MQ*tempMixing*(sman*xd0-8*xd00)+(2*xb0-sman*xc0+uman*xc0-sman*xc1-sman*xc2-tman*xc2+uman*xc2-MQ2*sman*xd0+4*MGl2*xd00+4*MQ2*xd00-2*tman*xd00-4*uman*xd00+MGl2*uman*xd1-MQ2*uman*xd1-tman*uman*xd12+uman**2*xd13+MGl2*sman*xd2-MQ2*sman*xd2-sman*tman*xd2+MGl2*uman*xd2-MQ2*uman*xd2-sman*tman*xd22-tman*uman*xd22
     &+sman*uman*xd23+uman**2*xd23+M2j**2*(xd1+2*xd11+3*xd12+xd13+xd2+xd22+xd23)-M2j*(xc0-2*xd00+MGl2*xd1-MQ2*xd1+sman*xd1+uman*xd1+2*uman*xd11+2*sman*xd12-tman*xd12+3*uman*xd12+2*uman*xd13+MGl2*xd2-MQ2*xd2+sman*xd2+uman*xd2+sman*xd22-tman*xd22+uman*xd22+sman*xd23+2*uman*xd23))*Krondelta(i,j)))/(8d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)-(tman*(xd23+xd3+xd33)+M2j*(2*xd13+xd23+xd3+xd33)-2*(MQ2*xd33+MGl2*(xd3+xd33)))*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd2+2*xd22)+(-(uman*xd2)+2*MGl2*xd22-uman*xd22+2*MQ2*(xd2+xd22)-uman*xd23+M2j*(2*xd12+xd2+xd22+xd23))*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd1+xd11+2*xd12+xd2+xd22)+(-2*xc0+4*xd00+4*MQ2*xd1+2*M2j*xd11+2*MGl2*xd11+2*MQ2*xd11+2*M2j*xd12+4*MGl2*xd12+4*MQ2*xd12+sman*xd12+2*tman*xd12+2*M2j*xd13-sman*xd13-2*uman*xd13+4*MQ2*xd2+sman*xd2+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+2*tman*xd22+2*M2j*xd23-sman*xd23-2*uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3)))+(xc0+xc1+xc2+MQ2*xd0+M2j*xd1+2*M2j*xd13+MGl2*xd2+MQ2*xd2+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-uman*xd23+M2j*xd3+2*MQ2*xd3-uman*xd3+M2j*xd33-uman*xd33)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0+xc2+2*MQ2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1+3*M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-uman*xd13+MGl2*xd2+MQ2*xd2+tman*xd2+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+2*tman*xd23-uman*xd23+2*M2j*xd3+4*MQ2*xd3-sman*xd3-2*uman*xd3+2*M2j*xd33-sman*xd33-2*uman*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(MGl*MQ*tempMixing*(xd0-4*xd23)-(xc0+xc1+xc2+MQ2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2-MGl2*xd2+MQ2*xd2+tman*xd2+M2j*xd22+tman*xd22+M2j*xd23-2*MGl2*xd23-2*MQ2*xd23+tman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-4*MGl*MQ*tempMixing*(xd12+xd2+xd22)+(2*xc1+xc2+2*xd00+3*M2j*xd12+2*MGl2*xd12+2*MQ2*xd12-uman*xd12+M2j*xd2+4*MQ2*xd2-uman*xd2+M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+2*tman*xd22-uman*xd22+2*M2j*xd23-sman*xd23-2*uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd13+xd23)+(xc0-xc2-4*xd00-M2j*xd1+MGl2*xd1-MQ2*xd1-2*M2j*xd11-3*M2j*xd12-tman*xd12-2*M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-tman*xd13+uman*xd13-M2j*xd2+MGl2*xd2-MQ2*xd2-M2j*xd22-tman*xd22-2*M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23-tman*xd23+uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd1+2*xd12+xd2+2*xd22)+(-2*xc0-2*xc1-xc2+2*xd00+2*MQ2*xd1+2*M2j*xd11+2*M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+M2j*xd13-uman*xd13+2*MQ2*xd2+sman*xd2+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2j*xd23-uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0/)      
      
      do n=1,17
      PA_A4(1,n)=PA_A4(1,n)+colors(1)*tempPA(n)
      enddo

!     Box 6.2 (28) DOES NOT AGREE WITH FORMCALC EXCEPT FOR Ag (1), (BECAUSE FEYNARTS AMPLITUDE IS WRONG)
      call LTDget(M2j,0d0,M2,0d0,tman,uman,MGl2,MQ2,MQ2,MGl2)
      call myCget(0d0,M2,uman,MQ2,MQ2,MGl2)    
      call myBget(M2,MQ2,MGl2)

      tempMixing=USf(i,2,type1,gen1)*USf(j,1,type1,gen1)+USf(i,1,type1,gen1)*USf(j,2,type1,gen1)
      colors=(/1d0/4d0,0d0,0d0/)
      tempPA=(/(GS**4*(MGl*MQ*tempMixing*(-(sman*xd0)+8*xd00)+(-2*xb0+sman*xc0-uman*xc0+sman*xc1+sman*xc2+tman*xc2-uman*xc2+MGl2*sman*xd0-4*MGl2*xd00-4*MQ2*xd00+2*tman*xd00+4*uman*xd00+MGl2*uman*xd1-MQ2*uman*xd1+tman*uman*xd12-uman**2*xd13+MGl2*sman*xd2-MQ2*sman*xd2+sman*tman*xd2+MGl2*uman*xd2-MQ2*uman*xd2+sman*tman*xd22+tman*uman*xd22
     &-sman*uman*xd23-uman**2*xd23-M2j**2*(xd1+2*xd11+3*xd12+xd13+xd2+xd22+xd23)+M2j*(xc0-2*xd00-MGl2*xd1+MQ2*xd1+sman*xd1+uman*xd1+2*uman*xd11+2*sman*xd12-tman*xd12+3*uman*xd12+2*uman*xd13-MGl2*xd2+MQ2*xd2+sman*xd2+uman*xd2+sman*xd22-tman*xd22+uman*xd22+sman*xd23+2*uman*xd23))*Krondelta(i,j)))/(8d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd3+2*xd33)-(tman*(xd23+xd3+xd33)+M2j*(2*xd13+xd23+xd3+xd33)-2*(MGl2*xd33+MQ2*(xd3+xd33)))*Krondelta(i,j)))/(4d0*Pi**2),
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd2+2*xd22)+(-(uman*xd2)+2*MQ2*xd22-uman*xd22+2*MGl2*(xd2+xd22)-uman*xd23+M2j*(2*xd12+xd2+xd22+xd23))*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd1+xd11+2*xd12+xd2+xd22)+(-2*xc0+4*xd00+4*MGl2*xd1+2*M2j*xd11+2*MGl2*xd11+2*MQ2*xd11+2*M2j*xd12+4*MGl2*xd12+4*MQ2*xd12+sman*xd12+2*tman*xd12+2*M2j*xd13-sman*xd13-2*uman*xd13+4*MGl2*xd2+sman*xd2+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+2*tman*xd22+2*M2j*xd23-sman*xd23-2*uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-(MGl*MQ*tempMixing*(xd0+2*(xd2+2*xd23+xd3)))+(xc0+xc1+xc2+MGl2*xd0+M2j*xd1+2*M2j*xd13+MGl2*xd2+MQ2*xd2+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23-uman*xd23+M2j*xd3+2*MGl2*xd3-uman*xd3+M2j*xd33-uman*xd33)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &-(GS**4*(-2*MGl*MQ*tempMixing*(xd0+xd1+2*xd13+xd2+2*xd23+2*xd3)+(xc0+xc2+2*MGl2*xd0+M2j*xd1+MGl2*xd1+MQ2*xd1+3*M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-uman*xd13+MGl2*xd2+MQ2*xd2+tman*xd2+M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+2*tman*xd23-uman*xd23+2*M2j*xd3+4*MGl2*xd3-sman*xd3-2*uman*xd3+2*M2j*xd33-sman*xd33-2*uman*xd33)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(MGl*MQ*tempMixing*(xd0-4*xd23)-(xc0+xc1+xc2+MGl2*xd0+M2j*xd1+2*M2j*xd12+M2j*xd2+MGl2*xd2-MQ2*xd2+tman*xd2+M2j*xd22+tman*xd22+M2j*xd23-2*MGl2*xd23-2*MQ2*xd23+tman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0,
     &(GS**4*(-4*MGl*MQ*tempMixing*(xd12+xd2+xd22)+(2*xc1+xc2+2*xd00+3*M2j*xd12+2*MGl2*xd12+2*MQ2*xd12-uman*xd12+M2j*xd2+4*MGl2*xd2-uman*xd2+M2j*xd22+2*MGl2*xd22+2*MQ2*xd22+2*tman*xd22-uman*xd22+2*M2j*xd23-sman*xd23-2*uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0,c0,c0,
     &-(GS**4*(-4*MGl*MQ*tempMixing*(xd13+xd23)+(xc0-xc2-4*xd00-M2j*xd1-MGl2*xd1+MQ2*xd1-2*M2j*xd11-3*M2j*xd12-tman*xd12-2*M2j*xd13+2*MGl2*xd13+2*MQ2*xd13-tman*xd13+uman*xd13-M2j*xd2-MGl2*xd2+MQ2*xd2-M2j*xd22-tman*xd22-2*M2j*xd23+2*MGl2*xd23+2*MQ2*xd23+sman*xd23-tman*xd23+uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),
     &(GS**4*(-2*MGl*MQ*tempMixing*(xd1+2*xd12+xd2+2*xd22)+(-2*xc0-2*xc1-xc2+2*xd00+2*MGl2*xd1+2*M2j*xd11+2*M2j*xd12+2*MGl2*xd12+2*MQ2*xd12+M2j*xd13-uman*xd13+2*MGl2*xd2+sman*xd2+2*MGl2*xd22+2*MQ2*xd22+sman*xd22+M2j*xd23-uman*xd23)*Krondelta(i,j)))/(4d0*Pi**2),c0/)      
      
      do n=1,17
      PA_A4(1,n)=PA_A4(1,n)+colors(1)*tempPA(n)
      enddo

      endif 

      end


      subroutine stst2gg_calc_M2Vertex4(stst2gg_M2Vertex4)
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double complex PA_A4(3,17)
      double precision M2,k1p1,k1p2,Krondelta,stst2gg_M2Vertex4(3)

      M2=MSf2(iSf1,itt1,igen1)

! Values for scalar products for the 2->2 process
      k1k2=sman/2d0
      p1p2=(sman-2*M2)/2d0
      k1p1=(M2-tman)/2d0
      k1p2=(M2-uman)/2d0

!Calculation of partial amplitudes and complex conjugation
      call stst2gg_calc_A4(PA_A4,iSf1,itt1,igen1,iSf2,itt2,igen2)
      PA_A4=dconjg(PA_A4)

! stst2gg_M2Vertex4 contains the singlet, symmetric octet and antisymmetric octet


! Singlet
      stst2gg_M2Vertex4(1) =   2*dble(-(gs**2*Nc*(-1 + Nc**2)*((-1 + FPSk1)*(-1 + FPSk2)*k1k2**4*(k1p1 + k1p2)*PA_A4(1,9) + 
     &       k1k2**3*(-((-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*(4*PA_A4(1,9) -PA_A4(1,11) - PA_A4(1,12))) + 
     &          k1p1*((-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*
     &              (4*PA_A4(1,2) + 4*PA_A4(1,3) + 4*PA_A4(1,6) + 4*PA_A4(1,9) - PA_A4(1,10) - PA_A4(1,11) - PA_A4(1,12) - PA_A4(1,15)) + 
     &             2*((-1 + FPSk1)*M2*PA_A4(1,10) + (-1 + FPSk1)*p1p2*PA_A4(1,11) + (-1 + FPSk2)*p1p2*PA_A4(1,12) + (-1 + FPSk2)*M2*PA_A4(1,15))) + 
     &          k1p2*((-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(4*PA_A4(1,9) - PA_A4(1,10) - PA_A4(1,15)) + 
     &             2*((-1 + FPSk1)*p1p2*PA_A4(1,10) + (-1 + FPSk1)*M2*PA_A4(1,11) + (-1 + FPSk2)*M2*PA_A4(1,12) + (-1 + FPSk2)*p1p2*PA_A4(1,15)))) + 
     &       4*FPSk1*FPSk2*(k1p1 + k1p2)*(k1p1**4*PA_A4(1,14) + k1p1*k1p2**3*(2*PA_A4(1,4) + 2*PA_A4(1,5) + PA_A4(1,14) - PA_A4(1,17)) + k1p2**4*PA_A4(1,17) + 
     &          k1p1**3*k1p2*(2*PA_A4(1,4) + 2*PA_A4(1,5) - PA_A4(1,14) + PA_A4(1,17)) + 3*k1p1**2*k1p2**2*(PA_A4(1,14) + PA_A4(1,17))) + 
     &       k1k2**2*(k1p1**3*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,9) - 2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A4(1,11) - 
     &             (-1 + FPSk2)*((-2 + 4*FPSk1)*PA_A4(1,12) - (-1 + FPSk1)*PA_A4(1,14))) + 
     &          k1p1**2*(4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(1,1) + 
     &             2*(p1p2*(-2*(-1 + FPSk1)*PA_A4(1,11) - 2*(-1 + FPSk2)*PA_A4(1,12) + (-2 + FPSk1 + FPSk2)*PA_A4(1,14)) + 
     &                M2*((-1 + FPSk1)*PA_A4(1,4) + (-1 + FPSk2)*PA_A4(1,5) - 2*(-1 + FPSk1)*PA_A4(1,10) - 2*(-1 + FPSk2)*PA_A4(1,15))) + 
     &             k1p2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,2) + 4*(-1 +FPSk1)*(-1 + FPSk2)*PA_A4(1,3) + PA_A4(1,4) - FPSk1*PA_A4(1,4) - FPSk2*PA_A4(1,4) + 
     &                FPSk1*FPSk2*PA_A4(1,4) + PA_A4(1,5) - FPSk1*PA_A4(1,5) - FPSk2*PA_A4(1,5) + FPSk1*FPSk2*PA_A4(1,5) + 4*PA_A4(1,6) - 4*FPSk1*PA_A4(1,6) - 
     &                4*FPSk2*PA_A4(1,6) + 4*FPSk1*FPSk2*PA_A4(1,6) - 4*PA_A4(1,7) + 4*FPSk1*PA_A4(1,7) + 4*FPSk2*PA_A4(1,7) - 4*FPSk1*FPSk2*PA_A4(1,7) - 
     &                4*PA_A4(1,8) + 4*FPSk1*PA_A4(1,8) + 4*FPSk2*PA_A4(1,8) - 4*FPSk1*FPSk2*PA_A4(1,8) - 2*PA_A4(1,10) + 2*FPSk1*PA_A4(1,10) + 8*FPSk2*PA_A4(1,10) - 
     &                8*FPSk1*FPSk2*PA_A4(1,10) - 4*PA_A4(1,11) + 4*FPSk1*PA_A4(1,11) + 4*FPSk2*PA_A4(1,11) - 4*FPSk1*FPSk2*PA_A4(1,11) - 4*PA_A4(1,12) + 
     &                4*FPSk1*PA_A4(1,12) + 4*FPSk2*PA_A4(1,12) - 4*FPSk1*FPSk2*PA_A4(1,12) - 4*PA_A4(1,13) + 4*FPSk1*PA_A4(1,13) + 4*FPSk2*PA_A4(1,13) - 
     &                4*FPSk1*FPSk2*PA_A4(1,13) + PA_A4(1,14) - FPSk1*PA_A4(1,14) - FPSk2*PA_A4(1,14) + FPSk1*FPSk2*PA_A4(1,14) - 2*PA_A4(1,15) + 
     &                8*FPSk1*PA_A4(1,15) + 2*FPSk2*PA_A4(1,15) - 8*FPSk1*FPSk2*PA_A4(1,15) - 4*PA_A4(1,16) + 4*FPSk1*PA_A4(1,16) + 4*FPSk2*PA_A4(1,16) - 
     &                4*FPSk1*FPSk2*PA_A4(1,16))) + k1p2*(4*(p1p2*(PA_A4(1,1) + M2*(PA_A4(1,4) + PA_A4(1,5))) + M2**2*PA_A4(1,14) + p1p2**2*PA_A4(1,17)) + 
     &             k1p2**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,9) - 2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A4(1,10) - 
     &                (-1 + FPSk2)*((-2 + 4*FPSk1)*PA_A4(1,15) - (-1 + FPSk1)*PA_A4(1,17))) + 
     &             2*k1p2*(2*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(1,1) + 
     &                M2*((-1 + FPSk2)*PA_A4(1,4) + (-1 + FPSk1)*PA_A4(1,5) - 2*(-1 + FPSk1)*PA_A4(1,11) - 2*(-1 + FPSk2)*PA_A4(1,12)) + 
     &                p1p2*(-2*(-1 + FPSk1)*PA_A4(1,10) - 2*(-1 + FPSk2)*PA_A4(1,15) + (-2 + FPSk1 + FPSk2)*PA_A4(1,17)))) + 
     &          k1p1*(k1p2**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,2) + 4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(1,3) + PA_A4(1,4) - FPSk1*PA_A4(1,4) - FPSk2*PA_A4(1,4) + 
     &                FPSk1*FPSk2*PA_A4(1,4) + PA_A4(1,5) - FPSk1*PA_A4(1,5) - FPSk2*PA_A4(1,5) + FPSk1*FPSk2*PA_A4(1,5) + 4*PA_A4(1,6) - 4*FPSk1*PA_A4(1,6) - 
     &                4*FPSk2*PA_A4(1,6) + 4*FPSk1*FPSk2*PA_A4(1,6) - 4*PA_A4(1,7) + 4*FPSk1*PA_A4(1,7) + 4*FPSk2*PA_A4(1,7) - 4*FPSk1*FPSk2*PA_A4(1,7) - 
     &                4*PA_A4(1,8) + 4*FPSk1*PA_A4(1,8) + 4*FPSk2*PA_A4(1,8) - 4*FPSk1*FPSk2*PA_A4(1,8) - 4*PA_A4(1,10) + 4*FPSk1*PA_A4(1,10) + 4*FPSk2*PA_A4(1,10) - 
     &                4*FPSk1*FPSk2*PA_A4(1,10) - 2*PA_A4(1,11) + 2*FPSk1*PA_A4(1,11) + 8*FPSk2*PA_A4(1,11) - 8*FPSk1*FPSk2*PA_A4(1,11) - 2*PA_A4(1,12) + 
     &                8*FPSk1*PA_A4(1,12) + 2*FPSk2*PA_A4(1,12) - 8*FPSk1*FPSk2*PA_A4(1,12) - 4*PA_A4(1,13) + 4*FPSk1*PA_A4(1,13) + 4*FPSk2*PA_A4(1,13) - 
     &                4*FPSk1*FPSk2*PA_A4(1,13) - 4*PA_A4(1,15) + 4*FPSk1*PA_A4(1,15) + 4*FPSk2*PA_A4(1,15) - 4*FPSk1*FPSk2*PA_A4(1,15) - 4*PA_A4(1,16) + 
     &                4*FPSk1*PA_A4(1,16) + 4*FPSk2*PA_A4(1,16) - 4*FPSk1*FPSk2*PA_A4(1,16) + PA_A4(1,17) - FPSk1*PA_A4(1,17) - FPSk2*PA_A4(1,17) + 
     &                FPSk1*FPSk2*PA_A4(1,17)) + 4*(p1p2*(PA_A4(1,1) + M2*(PA_A4(1,4) + PA_A4(1,5))) + p1p2**2*PA_A4(1,14) + M2**2*PA_A4(1,17)) + 
     &             2*k1p2*(4*(-2 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(1,1) + 
     &                M2*(-2*PA_A4(1,4) - 2*PA_A4(1,5) + 2*PA_A4(1,7) - 2*FPSk1*PA_A4(1,7) + 2*PA_A4(1,8) - 2*FPSk1*PA_A4(1,8) + 2*PA_A4(1,13) - 2*FPSk2*PA_A4(1,13) - 
     &                   2*PA_A4(1,14) + FPSk1*PA_A4(1,14) + FPSk2*PA_A4(1,14) + 2*PA_A4(1,16) - 2*FPSk2*PA_A4(1,16) - 2*PA_A4(1,17) + FPSk1*PA_A4(1,17) + 
     &                   FPSk2*PA_A4(1,17)) + p1p2*((-2 + FPSk1 + FPSk2)*PA_A4(1,4) + (-2 + FPSk1 + FPSk2)*PA_A4(1,5) - 
     &                   2*((-1 + FPSk1)*PA_A4(1,7) + (-1 + FPSk1)*PA_A4(1,8) - PA_A4(1,13) + FPSk2*PA_A4(1,13) + PA_A4(1,14) - PA_A4(1,16) + FPSk2*PA_A4(1,16) + 
     &                      PA_A4(1,17)))))) + 2*k1k2*(-(FPSk2*(k1p1**4*(2*PA_A4(1,11) - PA_A4(1,14)) + 
     &               k1p1**3*(2*PA_A4(1,1) + 2*(M2*PA_A4(1,5) + p1p2*PA_A4(1,14)) - 
     &                  k1p2*(3*PA_A4(1,4) + PA_A4(1,5) - 2*(PA_A4(1,7) + PA_A4(1,8) + 2*PA_A4(1,10) - PA_A4(1,14) - PA_A4(1,17)))) + 
     &               k1p1*k1p2**2*(2*PA_A4(1,1) + 2*(p1p2*(2*PA_A4(1,5) + PA_A4(1,14)) + M2*(PA_A4(1,5) + 2*PA_A4(1,14))) - 
     &                  k1p2*(PA_A4(1,4) + 3*PA_A4(1,5) - 2*(PA_A4(1,7) + PA_A4(1,8) + 2*PA_A4(1,11) - PA_A4(1,14) - PA_A4(1,17)))) + 
     &               k1p2**3*(2*PA_A4(1,1) + 2*M2*PA_A4(1,4) + 2*k1p2*PA_A4(1,10) - k1p2*PA_A4(1,17) + 2*p1p2*PA_A4(1,17)) + 
     &               k1p1**2*k1p2*(2*PA_A4(1,1) + k1p2*(-4*PA_A4(1,4) - 4*PA_A4(1,5) + 4*PA_A4(1,7) + 4*PA_A4(1,8) + 2*PA_A4(1,10) + 2*PA_A4(1,11) - 3*PA_A4(1,14) - 
     &                     3*PA_A4(1,17)) + 2*(M2*PA_A4(1,4) + 2*p1p2*PA_A4(1,4) + 2*M2*PA_A4(1,17) + p1p2*PA_A4(1,17))))) + 
     &          FPSk1*(k1p1**4*(-2*PA_A4(1,12) + 2*FPSk2*(PA_A4(1,11) + PA_A4(1,12) - PA_A4(1,14)) + PA_A4(1,14)) + 
     &             k1p2**3*(2*(-1 + FPSk2)*PA_A4(1,1) - 2*M2*PA_A4(1,5) + 2*FPSk2*k1p2*PA_A4(1,10) - 2*k1p2*PA_A4(1,15) + 2*FPSk2*k1p2*PA_A4(1,15) + 
     &                k1p2*PA_A4(1,17) - 2*FPSk2*k1p2*PA_A4(1,17) - 2*p1p2*PA_A4(1,17)) + 
     &             k1p1*k1p2**2*(2*(-1 + FPSk2)*PA_A4(1,1) - 2*(p1p2*(2*PA_A4(1,4) + PA_A4(1,14)) + M2*(PA_A4(1,4) + 2*PA_A4(1,14))) + 
     &                k1p2*((3 - 4*FPSk2)*PA_A4(1,4) + (1 - 4*FPSk2)*PA_A4(1,5) + 
     &                   2*(-2*PA_A4(1,12) - PA_A4(1,13) + PA_A4(1,14) - PA_A4(1,16) + 
     &                      FPSk2*(PA_A4(1,7) + PA_A4(1,8) + 2*PA_A4(1,11) + 2*PA_A4(1,12) + PA_A4(1,13) - PA_A4(1,14) + PA_A4(1,16) - PA_A4(1,17)) + PA_A4(1,17)))) + 
     &             k1p1**3*(2*(-1 + FPSk2)*PA_A4(1,1) - 2*(M2*PA_A4(1,4) +p1p2*PA_A4(1,14)) + 
     &                k1p2*((1 - 4*FPSk2)*PA_A4(1,4) + (3 - 4*FPSk2)*PA_A4(1,5) + 
     &                   2*(-PA_A4(1,13) + PA_A4(1,14) - 2*PA_A4(1,15) - PA_A4(1,16) + 
     &                      FPSk2*(PA_A4(1,7) + PA_A4(1,8) + 2*PA_A4(1,10) + PA_A4(1,13) - PA_A4(1,14) + 2*PA_A4(1,15) + PA_A4(1,16) - PA_A4(1,17)) + PA_A4(1,17)))) + 
     &             k1p1**2*k1p2*(2*(-1 + FPSk2)*PA_A4(1,1) + k1p2*
     &                 (-4*(-1 + FPSk2)*PA_A4(1,4) - 4*(-1 + FPSk2)*PA_A4(1,5) + 4*FPSk2*PA_A4(1,7) + 4*FPSk2*PA_A4(1,8) + 2*FPSk2*PA_A4(1,10) + 
     &                   2*FPSk2*PA_A4(1,11) - 2*PA_A4(1,12) + 2*FPSk2*PA_A4(1,12) - 4*PA_A4(1,13) + 4*FPSk2*PA_A4(1,13) + 3*PA_A4(1,14) - 6*FPSk2*PA_A4(1,14) - 
     &                   2*PA_A4(1,15) + 2*FPSk2*PA_A4(1,15) - 4*PA_A4(1,16) + 4*FPSk2*PA_A4(1,16) + 3*PA_A4(1,17) - 6*FPSk2*PA_A4(1,17)) - 
     &                2*(p1p2*(2*PA_A4(1,5) + PA_A4(1,17)) + M2*(PA_A4(1,5)+ 2*PA_A4(1,17))))))))/(12*k1k2**2*k1p1*k1p2))

      stst2gg_M2Vertex4(2) =  2*dble(-(CF*gs**2*(-4 + Nc**2)*((-1 + FPSk1)*(-1 + FPSk2)*k1k2**4*(k1p1 + k1p2)*PA_A4(2,9) + 
     &       k1k2**3*(-((-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*(4*PA_A4(2,9) - PA_A4(2,11) - PA_A4(2,12))) + 
     &          k1p1*((-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*
     &              (4*PA_A4(2,2) + 4*PA_A4(2,3) + 4*PA_A4(2,6) + 4*PA_A4(2,9) - PA_A4(2,10) - PA_A4(2,11) - PA_A4(2,12) - PA_A4(2,15)) + 
     &             2*((-1 + FPSk1)*M2*PA_A4(2,10) + (-1 + FPSk1)*p1p2*PA_A4(2,11) + (-1 + FPSk2)*p1p2*PA_A4(2,12) + (-1 + FPSk2)*M2*PA_A4(2,15))) + 
     &          k1p2*((-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(4*PA_A4(2,9) - PA_A4(2,10) - PA_A4(2,15)) + 
     &             2*((-1 + FPSk1)*p1p2*PA_A4(2,10) + (-1 + FPSk1)*M2*PA_A4(2,11) + (-1 + FPSk2)*M2*PA_A4(2,12) + (-1 + FPSk2)*p1p2*PA_A4(2,15)))) + 
     &       4*FPSk1*FPSk2*(k1p1 + k1p2)*(k1p1**4*PA_A4(2,14) + k1p1*k1p2**3*(2*PA_A4(2,4) + 2*PA_A4(2,5) + PA_A4(2,14) - PA_A4(2,17)) + k1p2**4*PA_A4(2,17) + 
     &          k1p1**3*k1p2*(2*PA_A4(2,4) + 2*PA_A4(2,5) - PA_A4(2,14) + PA_A4(2,17)) + 3*k1p1**2*k1p2**2*(PA_A4(2,14) + PA_A4(2,17))) + 
     &       k1k2**2*(k1p1**3*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,9) - 2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A4(2,11) - 
     &             (-1 + FPSk2)*((-2 + 4*FPSk1)*PA_A4(2,12) - (-1 + FPSk1)*PA_A4(2,14))) + 
     &          k1p1**2*(4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(2,1) + 
     &             2*(p1p2*(-2*(-1 + FPSk1)*PA_A4(2,11) - 2*(-1 + FPSk2)*PA_A4(2,12) + (-2 + FPSk1 + FPSk2)*PA_A4(2,14)) + 
     &                M2*((-1 + FPSk1)*PA_A4(2,4) + (-1 + FPSk2)*PA_A4(2,5) - 2*(-1 + FPSk1)*PA_A4(2,10) - 2*(-1 + FPSk2)*PA_A4(2,15))) + 
     &             k1p2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,2) + 4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,3) + PA_A4(2,4) - FPSk1*PA_A4(2,4) - FPSk2*PA_A4(2,4) + 
     &                FPSk1*FPSk2*PA_A4(2,4) + PA_A4(2,5) - FPSk1*PA_A4(2,5) - FPSk2*PA_A4(2,5) + FPSk1*FPSk2*PA_A4(2,5) + 4*PA_A4(2,6) - 4*FPSk1*PA_A4(2,6) - 
     &                4*FPSk2*PA_A4(2,6) + 4*FPSk1*FPSk2*PA_A4(2,6) - 4*PA_A4(2,7) + 4*FPSk1*PA_A4(2,7) + 4*FPSk2*PA_A4(2,7) - 4*FPSk1*FPSk2*PA_A4(2,7) - 
     &                4*PA_A4(2,8) + 4*FPSk1*PA_A4(2,8) + 4*FPSk2*PA_A4(2,8) - 4*FPSk1*FPSk2*PA_A4(2,8) - 2*PA_A4(2,10) + 2*FPSk1*PA_A4(2,10) + 8*FPSk2*PA_A4(2,10) - 
     &                8*FPSk1*FPSk2*PA_A4(2,10) - 4*PA_A4(2,11) + 4*FPSk1*PA_A4(2,11) + 4*FPSk2*PA_A4(2,11) - 4*FPSk1*FPSk2*PA_A4(2,11) - 4*PA_A4(2,12) + 
     &                4*FPSk1*PA_A4(2,12) + 4*FPSk2*PA_A4(2,12) - 4*FPSk1*FPSk2*PA_A4(2,12) - 4*PA_A4(2,13) + 4*FPSk1*PA_A4(2,13) + 4*FPSk2*PA_A4(2,13) - 
     &                4*FPSk1*FPSk2*PA_A4(2,13) + PA_A4(2,14) - FPSk1*PA_A4(2,14) - FPSk2*PA_A4(2,14) + FPSk1*FPSk2*PA_A4(2,14) - 2*PA_A4(2,15) + 
     &                8*FPSk1*PA_A4(2,15) + 2*FPSk2*PA_A4(2,15) - 8*FPSk1*FPSk2*PA_A4(2,15) - 4*PA_A4(2,16) + 4*FPSk1*PA_A4(2,16) + 4*FPSk2*PA_A4(2,16) - 
     &                4*FPSk1*FPSk2*PA_A4(2,16))) + k1p2*(4*(p1p2*(PA_A4(2,1) + M2*(PA_A4(2,4) + PA_A4(2,5))) + M2**2*PA_A4(2,14) + p1p2**2*PA_A4(2,17)) + 
     &             k1p2**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,9) - 2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A4(2,10) - 
     &                (-1 + FPSk2)*((-2 + 4*FPSk1)*PA_A4(2,15) - (-1 + FPSk1)*PA_A4(2,17))) + 
     &             2*k1p2*(2*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(2,1) + 
     &                M2*((-1 + FPSk2)*PA_A4(2,4) + (-1 + FPSk1)*PA_A4(2,5) - 2*(-1 + FPSk1)*PA_A4(2,11) - 2*(-1 + FPSk2)*PA_A4(2,12)) + 
     &                p1p2*(-2*(-1 + FPSk1)*PA_A4(2,10) - 2*(-1 + FPSk2)*PA_A4(2,15) + (-2 + FPSk1 + FPSk2)*PA_A4(2,17)))) + 
     &          k1p1*(k1p2**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,2) + 4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(2,3) + PA_A4(2,4) - FPSk1*PA_A4(2,4) - FPSk2*PA_A4(2,4) + 
     &                FPSk1*FPSk2*PA_A4(2,4) + PA_A4(2,5) - FPSk1*PA_A4(2,5) - FPSk2*PA_A4(2,5) + FPSk1*FPSk2*PA_A4(2,5) + 4*PA_A4(2,6) - 4*FPSk1*PA_A4(2,6) - 
     &                4*FPSk2*PA_A4(2,6) + 4*FPSk1*FPSk2*PA_A4(2,6) - 4*PA_A4(2,7) + 4*FPSk1*PA_A4(2,7) + 4*FPSk2*PA_A4(2,7) - 4*FPSk1*FPSk2*PA_A4(2,7) - 
     &                4*PA_A4(2,8) + 4*FPSk1*PA_A4(2,8) + 4*FPSk2*PA_A4(2,8) - 4*FPSk1*FPSk2*PA_A4(2,8) - 4*PA_A4(2,10) + 4*FPSk1*PA_A4(2,10) + 4*FPSk2*PA_A4(2,10) - 
     &                4*FPSk1*FPSk2*PA_A4(2,10) - 2*PA_A4(2,11) + 2*FPSk1*PA_A4(2,11) + 8*FPSk2*PA_A4(2,11) - 8*FPSk1*FPSk2*PA_A4(2,11) - 2*PA_A4(2,12) + 
     &                8*FPSk1*PA_A4(2,12) + 2*FPSk2*PA_A4(2,12) - 8*FPSk1*FPSk2*PA_A4(2,12) - 4*PA_A4(2,13) + 4*FPSk1*PA_A4(2,13) + 4*FPSk2*PA_A4(2,13) - 
     &                4*FPSk1*FPSk2*PA_A4(2,13) - 4*PA_A4(2,15) + 4*FPSk1*PA_A4(2,15) + 4*FPSk2*PA_A4(2,15) - 4*FPSk1*FPSk2*PA_A4(2,15) - 4*PA_A4(2,16) + 
     &                4*FPSk1*PA_A4(2,16) + 4*FPSk2*PA_A4(2,16) - 4*FPSk1*FPSk2*PA_A4(2,16) + PA_A4(2,17) - FPSk1*PA_A4(2,17) - FPSk2*PA_A4(2,17) + 
     &                FPSk1*FPSk2*PA_A4(2,17)) + 4*(p1p2*(PA_A4(2,1) + M2*(PA_A4(2,4) + PA_A4(2,5))) + p1p2**2*PA_A4(2,14) + M2**2*PA_A4(2,17)) + 
     &             2*k1p2*(4*(-2 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(2,1) + 
     &                M2*(-2*PA_A4(2,4) - 2*PA_A4(2,5) + 2*PA_A4(2,7) - 2*FPSk1*PA_A4(2,7) + 2*PA_A4(2,8) - 2*FPSk1*PA_A4(2,8) + 2*PA_A4(2,13) - 2*FPSk2*PA_A4(2,13) - 
     &                   2*PA_A4(2,14) + FPSk1*PA_A4(2,14) + FPSk2*PA_A4(2,14) + 2*PA_A4(2,16) - 2*FPSk2*PA_A4(2,16) - 2*PA_A4(2,17) + FPSk1*PA_A4(2,17) + 
     &                   FPSk2*PA_A4(2,17)) + p1p2*((-2 + FPSk1 + FPSk2)*PA_A4(2,4) + (-2 + FPSk1 + FPSk2)*PA_A4(2,5) - 
     &                   2*((-1 + FPSk1)*PA_A4(2,7) + (-1 + FPSk1)*PA_A4(2,8) - PA_A4(2,13) + FPSk2*PA_A4(2,13) + PA_A4(2,14) - PA_A4(2,16) + FPSk2*PA_A4(2,16) + 
     &                      PA_A4(2,17)))))) + 2*k1k2*(-(FPSk2*(k1p1**4*(2*PA_A4(2,11) - PA_A4(2,14)) + 
     &               k1p1**3*(2*PA_A4(2,1) + 2*(M2*PA_A4(2,5) + p1p2*PA_A4(2,14)) - 
     &                  k1p2*(3*PA_A4(2,4) + PA_A4(2,5) - 2*(PA_A4(2,7) + PA_A4(2,8) + 2*PA_A4(2,10) - PA_A4(2,14) - PA_A4(2,17)))) + 
     &               k1p1*k1p2**2*(2*PA_A4(2,1) + 2*(p1p2*(2*PA_A4(2,5) + PA_A4(2,14)) + M2*(PA_A4(2,5) + 2*PA_A4(2,14))) - 
     &                  k1p2*(PA_A4(2,4) + 3*PA_A4(2,5) - 2*(PA_A4(2,7) + PA_A4(2,8) + 2*PA_A4(2,11) - PA_A4(2,14) - PA_A4(2,17)))) + 
     &               k1p2**3*(2*PA_A4(2,1) + 2*M2*PA_A4(2,4) + 2*k1p2*PA_A4(2,10) - k1p2*PA_A4(2,17) + 2*p1p2*PA_A4(2,17)) + 
     &               k1p1**2*k1p2*(2*PA_A4(2,1) + k1p2*(-4*PA_A4(2,4) - 4*PA_A4(2,5) + 4*PA_A4(2,7) + 4*PA_A4(2,8) + 2*PA_A4(2,10) + 2*PA_A4(2,11) - 3*PA_A4(2,14) - 
     &                     3*PA_A4(2,17)) + 2*(M2*PA_A4(2,4) + 2*p1p2*PA_A4(2,4) + 2*M2*PA_A4(2,17) + p1p2*PA_A4(2,17))))) + 
     &          FPSk1*(k1p1**4*(-2*PA_A4(2,12) + 2*FPSk2*(PA_A4(2,11) + PA_A4(2,12) - PA_A4(2,14)) + PA_A4(2,14)) + 
     &             k1p2**3*(2*(-1 + FPSk2)*PA_A4(2,1) - 2*M2*PA_A4(2,5) + 2*FPSk2*k1p2*PA_A4(2,10) - 2*k1p2*PA_A4(2,15) + 2*FPSk2*k1p2*PA_A4(2,15) + 
     &                k1p2*PA_A4(2,17) - 2*FPSk2*k1p2*PA_A4(2,17) - 2*p1p2*PA_A4(2,17)) + 
     &             k1p1*k1p2**2*(2*(-1 + FPSk2)*PA_A4(2,1) - 2*(p1p2*(2*PA_A4(2,4) + PA_A4(2,14)) + M2*(PA_A4(2,4) + 2*PA_A4(2,14))) + 
     &                k1p2*((3 - 4*FPSk2)*PA_A4(2,4) + (1 - 4*FPSk2)*PA_A4(2,5) + 
     &                   2*(-2*PA_A4(2,12) - PA_A4(2,13) + PA_A4(2,14) - PA_A4(2,16) + 
     &                      FPSk2*(PA_A4(2,7) + PA_A4(2,8) + 2*PA_A4(2,11) + 2*PA_A4(2,12) + PA_A4(2,13) - PA_A4(2,14) + PA_A4(2,16) - PA_A4(2,17)) + PA_A4(2,17)))) + 
     &             k1p1**3*(2*(-1 + FPSk2)*PA_A4(2,1) - 2*(M2*PA_A4(2,4) + p1p2*PA_A4(2,14)) + 
     &                k1p2*((1 - 4*FPSk2)*PA_A4(2,4) + (3 - 4*FPSk2)*PA_A4(2,5) + 
     &                   2*(-PA_A4(2,13) + PA_A4(2,14) - 2*PA_A4(2,15) - PA_A4(2,16) + 
     &                      FPSk2*(PA_A4(2,7) + PA_A4(2,8) + 2*PA_A4(2,10) + PA_A4(2,13) - PA_A4(2,14) + 2*PA_A4(2,15) + PA_A4(2,16) - PA_A4(2,17)) + PA_A4(2,17)))) + 
     &             k1p1**2*k1p2*(2*(-1 + FPSk2)*PA_A4(2,1) + k1p2*
     &                 (-4*(-1 + FPSk2)*PA_A4(2,4) - 4*(-1 + FPSk2)*PA_A4(2,5) + 4*FPSk2*PA_A4(2,7) + 4*FPSk2*PA_A4(2,8) + 2*FPSk2*PA_A4(2,10) + 
     &                   2*FPSk2*PA_A4(2,11) - 2*PA_A4(2,12) + 2*FPSk2*PA_A4(2,12) - 4*PA_A4(2,13) + 4*FPSk2*PA_A4(2,13) + 3*PA_A4(2,14) - 6*FPSk2*PA_A4(2,14) - 
     &                   2*PA_A4(2,15) + 2*FPSk2*PA_A4(2,15) - 4*PA_A4(2,16) + 4*FPSk2*PA_A4(2,16) + 3*PA_A4(2,17) - 6*FPSk2*PA_A4(2,17)) - 
     &                2*(p1p2*(2*PA_A4(2,5) + PA_A4(2,17)) + M2*(PA_A4(2,5) + 2*PA_A4(2,17))))))))/(4*k1k2**2*k1p1*k1p2))


      stst2gg_M2Vertex4(3) =  2*dble((CF*CV*gs**2*Nc*(-((-1 + FPSk1)*(-1 + FPSk2)*k1k2**4*(k1p1 - k1p2)*sman*PA_A4(3,9)) + 
     &      k1k2**3*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*(4*k1p2*(PA_A4(3,2) + PA_A4(3,3) - 2*PA_A4(3,6) + 4*PA_A4(3,9)) + 
     &            sman*(4*PA_A4(3,9) - PA_A4(3,11) - PA_A4(3,12))) + 
     &         k1p2*sman*((-1 + FPSk1)*(-1 + FPSk2)*PA_A4(3,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(4*PA_A4(3,9) - PA_A4(3,10) - PA_A4(3,15)) + 
     &            2*((-1 + FPSk1)*p1p2*PA_A4(3,10) + (-1 + FPSk1)*M2*PA_A4(3,11) + (-1 + FPSk2)*M2*PA_A4(3,12) + (-1 + FPSk2)*p1p2*PA_A4(3,15))) + 
     &         k1p1*(sman*((-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,1) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(PA_A4(3,10) - PA_A4(3,11) - PA_A4(3,12) + PA_A4(3,15)) - 
     &               2*((-1 + FPSk1)*M2*PA_A4(3,10) + (-1 + FPSk1)*p1p2*PA_A4(3,11) + (-1 + FPSk2)*p1p2*PA_A4(3,12) + (-1 + FPSk2)*M2*PA_A4(3,15))) - 
     &            4*k1p2*((-1 + FPSk1)*(-1 + FPSk2)*k1p2*(PA_A4(3,2) + PA_A4(3,3) - 2*PA_A4(3,6) + 4*PA_A4(3,9)) + 
     &               (M2 - p1p2)*(2*(-1 + FPSk1)*PA_A4(3,7) - 2*(-1 + FPSk1)*PA_A4(3,8) - PA_A4(3,10) + FPSk1*PA_A4(3,10) + PA_A4(3,11) - FPSk1*PA_A4(3,11) + 
     &                  PA_A4(3,12) - FPSk2*PA_A4(3,12) + 2*PA_A4(3,13) - 2*FPSk2*PA_A4(3,13) - PA_A4(3,15) + FPSk2*PA_A4(3,15) - 2*PA_A4(3,16) + 2*FPSk2*PA_A4(3,16))))
     &         ) - 4*FPSk1*FPSk2*(k1p1 - k1p2)*sman*(k1p1**4*PA_A4(3,14) + k1p1**3*k1p2*(2*PA_A4(3,4) + 2*PA_A4(3,5) + PA_A4(3,14) - PA_A4(3,17)) + 
     &         k1p2**4*PA_A4(3,17) + k1p1*k1p2**3*(2*PA_A4(3,4) + 2*PA_A4(3,5) - PA_A4(3,14) + PA_A4(3,17)) + 3*k1p1**2*k1p2**2*(PA_A4(3,14) + PA_A4(3,17))) - 
     &      2*k1k2*(FPSk1*(k1p2**3*sman*(-2*(-1 + FPSk2)*PA_A4(3,1) + 2*M2*PA_A4(3,5) - 2*FPSk2*k1p2*PA_A4(3,10) + 2*k1p2*PA_A4(3,15) - 2*FPSk2*k1p2*PA_A4(3,15) - 
     &               k1p2*PA_A4(3,17) + 2*FPSk2*k1p2*PA_A4(3,17) + 2*p1p2*PA_A4(3,17)) + 
     &            k1p1**3*(sman*(2*(-1 + FPSk2)*PA_A4(3,1) - 2*(M2*PA_A4(3,4) + p1p2*PA_A4(3,14)) + 
     &                  k1p2*((1 - 2*FPSk2)*PA_A4(3,4) + (1 - 2*FPSk2)*PA_A4(3,5) + 
     &                     2*(-PA_A4(3,13) - 2*PA_A4(3,15) + FPSk2*(-PA_A4(3,7) + PA_A4(3,8) + 2*PA_A4(3,10) + PA_A4(3,13) + 2*PA_A4(3,15) - PA_A4(3,16)) + PA_A4(3,16)))
     &                  ) + 2*k1p2**2*(-((-1 + FPSk2)*PA_A4(3,4)) - (-4 + FPSk2)*PA_A4(3,5) - 2*PA_A4(3,14) + 2*FPSk2*PA_A4(3,14) + PA_A4(3,17) - 4*FPSk2*PA_A4(3,17))
     &               ) + k1p1**4*(sman*(-2*PA_A4(3,12) + 2*FPSk2*(PA_A4(3,11) + PA_A4(3,12) - PA_A4(3,14)) + PA_A4(3,14)) - 
     &               2*k1p2*((-2 + FPSk2)*PA_A4(3,4) + (1 + FPSk2)*PA_A4(3,5) - 3*PA_A4(3,14) + 4*FPSk2*PA_A4(3,14) - 2*FPSk2*PA_A4(3,17))) + 
     &            k1p1*k1p2**2*(sman*(2*(-1 + FPSk2)*PA_A4(3,1) - 2*M2*PA_A4(3,4) + 4*p1p2*PA_A4(3,4) + 4*M2*PA_A4(3,14) - 2*p1p2*PA_A4(3,14) + 
     &                  k1p2*((-1 + 2*FPSk2)*PA_A4(3,4) + (-1 + 2*FPSk2)*PA_A4(3,5) + 4*PA_A4(3,12) - 2*PA_A4(3,13) + 2*PA_A4(3,16) - 
     &                     2*FPSk2*(PA_A4(3,7) - PA_A4(3,8) + 2*PA_A4(3,11) + 2*PA_A4(3,12) - PA_A4(3,13) + PA_A4(3,16)))) + 
     &               2*k1p2**2*((1 + FPSk2)*PA_A4(3,4) + (-2 + FPSk2)*PA_A4(3,5) - 2*FPSk2*PA_A4(3,14) - 3*PA_A4(3,17) + 4*FPSk2*PA_A4(3,17))) + 
     &            k1p1**2*k1p2*(2*k1p2**2*((-4 + FPSk2)*PA_A4(3,4) + (-1 + FPSk2)*PA_A4(3,5) - PA_A4(3,14) + 4*FPSk2*PA_A4(3,14) + 2*PA_A4(3,17) - 
     &                  2*FPSk2*PA_A4(3,17)) + sman*(-2*(-1 + FPSk2)*PA_A4(3,1) - 4*p1p2*PA_A4(3,5) + 4*FPSk2*k1p2*PA_A4(3,7) - 4*FPSk2*k1p2*PA_A4(3,8) - 
     &                  2*FPSk2*k1p2*PA_A4(3,10) + 2*FPSk2*k1p2*PA_A4(3,11) - 2*k1p2*PA_A4(3,12) + 2*FPSk2*k1p2*PA_A4(3,12) + 4*k1p2*PA_A4(3,13) - 
     &                  4*FPSk2*k1p2*PA_A4(3,13) - k1p2*PA_A4(3,14) + 2*FPSk2*k1p2*PA_A4(3,14) + 2*k1p2*PA_A4(3,15) - 2*FPSk2*k1p2*PA_A4(3,15) - 4*k1p2*PA_A4(3,16) + 
     &                  4*FPSk2*k1p2*PA_A4(3,16) + 2*M2*(PA_A4(3,5) - 2*PA_A4(3,17)) + k1p2*PA_A4(3,17) - 2*FPSk2*k1p2*PA_A4(3,17) + 2*p1p2*PA_A4(3,17)))) + 
     &         FPSk2*(k1p1**4*(sman*(-2*PA_A4(3,11) + PA_A4(3,14)) + k1p2*(-2*PA_A4(3,4) + 4*PA_A4(3,5) + 6*PA_A4(3,14))) + 
     &            k1p2**3*sman*(2*PA_A4(3,1) + 2*M2*PA_A4(3,4) + 2*k1p2*PA_A4(3,10) - k1p2*PA_A4(3,17) + 2*p1p2*PA_A4(3,17)) + 
     &            k1p1**3*(sman*(-2*PA_A4(3,1) + k1p2*(PA_A4(3,4) + PA_A4(3,5) + 2*PA_A4(3,7) - 2*PA_A4(3,8) - 4*PA_A4(3,10)) - 2*(M2*PA_A4(3,5) + p1p2*PA_A4(3,14))) + 
     &               2*k1p2**2*(4*PA_A4(3,4) + PA_A4(3,5) - 2*PA_A4(3,14) + PA_A4(3,17))) - 
     &            k1p1*k1p2**2*(k1p2*sman*(PA_A4(3,4) + PA_A4(3,5) - 2*PA_A4(3,7) + 2*PA_A4(3,8) - 4*PA_A4(3,11)) + 
     &               2*sman*(PA_A4(3,1) + M2*(PA_A4(3,5) - 2*PA_A4(3,14)) + p1p2*(-2*PA_A4(3,5) + PA_A4(3,14))) + k1p2**2*(4*PA_A4(3,4) - 2*PA_A4(3,5) + 6*PA_A4(3,17)))
     &              + k1p1**2*k1p2*(-2*k1p2**2*(PA_A4(3,4) + 4*PA_A4(3,5) + PA_A4(3,14) - 2*PA_A4(3,17)) + 
     &               k1p2*sman*(-4*PA_A4(3,7) + 4*PA_A4(3,8) + 2*PA_A4(3,10) - 2*PA_A4(3,11) - PA_A4(3,14) + PA_A4(3,17)) + 
     &               2*sman*(PA_A4(3,1) + M2*(PA_A4(3,4) - 2*PA_A4(3,17)) + p1p2*(-2*PA_A4(3,4) + PA_A4(3,17)))))) + 
     &      k1k2**2*(k1p1**3*(sman*(4*(-1 + FPSk1)*PA_A4(3,9) - 2*(-1 + FPSk1)*PA_A4(3,11) + 2*PA_A4(3,12) - 4*FPSk1*PA_A4(3,12) - PA_A4(3,14) + FPSk1*PA_A4(3,14)) + 
     &            4*k1p2*((-1 + FPSk1)*PA_A4(3,8) - 2*(-1 + FPSk1)*PA_A4(3,10) + 3*PA_A4(3,11) - 3*FPSk1*PA_A4(3,11) + 3*PA_A4(3,12) - 4*FPSk1*PA_A4(3,12) - 
     &               PA_A4(3,13) - FPSk1*PA_A4(3,13) + 2*PA_A4(3,15) - FPSk1*PA_A4(3,15) + 2*FPSk1*PA_A4(3,16)) + 
     &            FPSk2*(sman*(-4*(-1 + FPSk1)*PA_A4(3,9) + 4*(-1 + FPSk1)*PA_A4(3,11) - 2*PA_A4(3,12) + 4*FPSk1*PA_A4(3,12) + PA_A4(3,14) - FPSk1*PA_A4(3,14)) - 
     &               4*k1p2*(2*(-1 + FPSk1)*PA_A4(3,7) - (-1 + FPSk1)*PA_A4(3,8) + PA_A4(3,10) - FPSk1*PA_A4(3,10) + 4*PA_A4(3,11) - 4*FPSk1*PA_A4(3,11) + 
     &                  3*PA_A4(3,12) - 4*FPSk1*PA_A4(3,12) - PA_A4(3,13) - FPSk1*PA_A4(3,13) + 2*PA_A4(3,15) - FPSk1*PA_A4(3,15) + 2*FPSk1*PA_A4(3,16)))) + 
     &         k1p2*sman*(4*(p1p2*(PA_A4(3,1) + M2*(PA_A4(3,4) + PA_A4(3,5))) + M2**2*PA_A4(3,14) + p1p2**2*PA_A4(3,17)) + 
     &            k1p2**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(3,9) - 2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A4(3,10) - 
     &               (-1 + FPSk2)*((-2 + 4*FPSk1)*PA_A4(3,15) - (-1 + FPSk1)*PA_A4(3,17))) + 
     &            2*k1p2*(2*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,1) + 
     &               M2*((-1 + FPSk2)*PA_A4(3,4) + (-1 + FPSk1)*PA_A4(3,5) - 2*(-1 + FPSk1)*PA_A4(3,11) - 2*(-1 + FPSk2)*PA_A4(3,12)) + 
     &               p1p2*(-2*(-1 + FPSk1)*PA_A4(3,10) - 2*(-1 + FPSk2)*PA_A4(3,15) + (-2 + FPSk1 + FPSk2)*PA_A4(3,17)))) - 
     &         k1p1*(4*k1p2**3*((-1 + FPSk1)*(1 + FPSk2)*PA_A4(3,7) + 3*PA_A4(3,10) - 3*FPSk1*PA_A4(3,10) + 2*PA_A4(3,11) - 2*FPSk1*PA_A4(3,11) + 2*PA_A4(3,12) - 
     &               FPSk1*PA_A4(3,12) + 2*FPSk1*PA_A4(3,13) + 3*PA_A4(3,15) - 4*FPSk1*PA_A4(3,15) - PA_A4(3,16) - FPSk1*PA_A4(3,16) + 
     &               FPSk2*(-2*(-1 + FPSk1)*PA_A4(3,8) + 4*(-1 + FPSk1)*PA_A4(3,10) - PA_A4(3,11) + FPSk1*PA_A4(3,11) - 2*PA_A4(3,12) + FPSk1*PA_A4(3,12) - 
     &                  2*FPSk1*PA_A4(3,13) - 3*PA_A4(3,15) + 4*FPSk1*PA_A4(3,15) + PA_A4(3,16) + FPSk1*PA_A4(3,16))) - 
     &            2*k1p2*sman*(FPSk1*p1p2*(PA_A4(3,4) - PA_A4(3,5) - 2*PA_A4(3,7) + 2*PA_A4(3,8)) + 
     &               FPSk2*p1p2*(-PA_A4(3,4) + PA_A4(3,5) + 2*PA_A4(3,13) - 2*PA_A4(3,16)) + 2*p1p2*(PA_A4(3,7) - PA_A4(3,8) - PA_A4(3,13) + PA_A4(3,16)) + 
     &               FPSk1*M2*(2*PA_A4(3,7) - 2*PA_A4(3,8) + PA_A4(3,14) - PA_A4(3,17)) - 
     &               2*M2*(PA_A4(3,7) - PA_A4(3,8) - PA_A4(3,13) + PA_A4(3,14) + PA_A4(3,16) - PA_A4(3,17)) + 
     &               FPSk2*M2*(-2*PA_A4(3,13) + PA_A4(3,14) + 2*PA_A4(3,16) - PA_A4(3,17))) + 
     &            4*sman*(p1p2*(PA_A4(3,1) + M2*(PA_A4(3,4) + PA_A4(3,5))) + p1p2**2*PA_A4(3,14) + M2**2*PA_A4(3,17)) + 
     &            k1p2**2*(8*(3 + FPSk1*(-1 + FPSk2) - FPSk2)*PA_A4(3,1) + 
     &               sman*(4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,2) + 4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,3) - PA_A4(3,4) + FPSk1*PA_A4(3,4) + 
     &                  FPSk2*PA_A4(3,4) - FPSk1*FPSk2*PA_A4(3,4) - PA_A4(3,5) + FPSk1*PA_A4(3,5) + FPSk2*PA_A4(3,5) - FPSk1*FPSk2*PA_A4(3,5) + 4*PA_A4(3,6) - 
     &                  4*FPSk1*PA_A4(3,6) - 4*FPSk2*PA_A4(3,6) + 4*FPSk1*FPSk2*PA_A4(3,6) + 2*PA_A4(3,11) - 2*FPSk1*PA_A4(3,11) - 4*FPSk2*PA_A4(3,11) + 
     &                  4*FPSk1*FPSk2*PA_A4(3,11) + 2*PA_A4(3,12) - 4*FPSk1*PA_A4(3,12) - 2*FPSk2*PA_A4(3,12) + 4*FPSk1*FPSk2*PA_A4(3,12) + PA_A4(3,17) - 
     &                  FPSk1*PA_A4(3,17) - FPSk2*PA_A4(3,17) + FPSk1*FPSk2*PA_A4(3,17)) + 
     &               4*(M2*((1 + 2*FPSk1 - FPSk2)*PA_A4(3,4) + (1 - FPSk1 + 2*FPSk2)*PA_A4(3,5) - (-2 + FPSk1 + FPSk2)*(2*PA_A4(3,14) - PA_A4(3,17))) + 
     &                  p1p2*((1 - 2*FPSk1 + FPSk2)*PA_A4(3,4) + (1 + FPSk1 - 2*FPSk2)*PA_A4(3,5) - 2*PA_A4(3,14) + 2*FPSk1*PA_A4(3,14) + 2*FPSk2*PA_A4(3,14) + 
     &                     4*PA_A4(3,17) - FPSk1*PA_A4(3,17) - FPSk2*PA_A4(3,17))))) + 
     &         k1p1**2*(2*sman*(2*(-1 + FPSk1)*(-1 + FPSk2)*PA_A4(3,1) + 
     &               p1p2*(2*(-1 + FPSk1)*PA_A4(3,11) + 2*(-1 + FPSk2)*PA_A4(3,12) - (-2 + FPSk1 + FPSk2)*PA_A4(3,14)) + 
     &               M2*(-((-1 + FPSk1)*PA_A4(3,4)) + PA_A4(3,5) - FPSk2*PA_A4(3,5) - 2*PA_A4(3,10) + 2*FPSk1*PA_A4(3,10) - 2*PA_A4(3,15) + 2*FPSk2*PA_A4(3,15))) + 
     &            4*k1p2**2*((-1 + FPSk1)*(1 + 3*FPSk2)*PA_A4(3,7) - (-1 + FPSk1)*(1 + 3*FPSk2)*PA_A4(3,8) + PA_A4(3,10) - FPSk1*PA_A4(3,10) - 3*FPSk2*PA_A4(3,10) + 
     &               3*FPSk1*FPSk2*PA_A4(3,10) - PA_A4(3,11) + FPSk1*PA_A4(3,11) + 3*FPSk2*PA_A4(3,11) - 3*FPSk1*FPSk2*PA_A4(3,11) - PA_A4(3,12) + 
     &               3*FPSk1*PA_A4(3,12) + FPSk2*PA_A4(3,12) - 3*FPSk1*FPSk2*PA_A4(3,12) + PA_A4(3,13) + 3*FPSk1*PA_A4(3,13) - FPSk2*PA_A4(3,13) - 
     &               3*FPSk1*FPSk2*PA_A4(3,13) + PA_A4(3,15) - 3*FPSk1*PA_A4(3,15) - FPSk2*PA_A4(3,15) + 3*FPSk1*FPSk2*PA_A4(3,15) - PA_A4(3,16) - 
     &               3*FPSk1*PA_A4(3,16) + FPSk2*PA_A4(3,16) + 3*FPSk1*FPSk2*PA_A4(3,16)) + 
     &            k1p2*(8*(3 + FPSk1*(-1 + FPSk2) - FPSk2)*PA_A4(3,1) + 
     &               sman*(4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,2) + 4*(-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*PA_A4(3,3) - PA_A4(3,4) + FPSk1*PA_A4(3,4) + 
     &                  FPSk2*PA_A4(3,4) - FPSk1*FPSk2*PA_A4(3,4) - PA_A4(3,5) + FPSk1*PA_A4(3,5) + FPSk2*PA_A4(3,5) - FPSk1*FPSk2*PA_A4(3,5) + 4*PA_A4(3,6) - 
     &                  4*FPSk1*PA_A4(3,6) - 4*FPSk2*PA_A4(3,6) + 4*FPSk1*FPSk2*PA_A4(3,6) + 2*PA_A4(3,10) - 2*FPSk1*PA_A4(3,10) - 4*FPSk2*PA_A4(3,10) + 
     &                  4*FPSk1*FPSk2*PA_A4(3,10) + PA_A4(3,14) - FPSk1*PA_A4(3,14) - FPSk2*PA_A4(3,14) + FPSk1*FPSk2*PA_A4(3,14) + 2*PA_A4(3,15) - 
     &                  4*FPSk1*PA_A4(3,15) - 2*FPSk2*PA_A4(3,15) + 4*FPSk1*FPSk2*PA_A4(3,15)) + 
     &               4*(M2*((1 - FPSk1 + 2*FPSk2)*PA_A4(3,4) + (1 + 2*FPSk1 - FPSk2)*PA_A4(3,5) + (-2 + FPSk1 + FPSk2)*(PA_A4(3,14) - 2*PA_A4(3,17))) + 
     &                  p1p2*((1 + FPSk1 - 2*FPSk2)*PA_A4(3,4) + (1 - 2*FPSk1 + FPSk2)*PA_A4(3,5) + 4*PA_A4(3,14) - FPSk1*PA_A4(3,14) - FPSk2*PA_A4(3,14) - 
     &                     2*PA_A4(3,17) + 2*FPSk1*PA_A4(3,17) + 2*FPSk2*PA_A4(3,17))))))))/(4*k1k2**2*k1p1*k1p2*sman))

      end