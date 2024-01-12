!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpFF(channel_,f3_type,f3_gen,f4_type,f4_gen,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,channel_,f3_type,f3_gen,f4_type,f4_gen
      double precision FF_TX(6),FF_UX(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox

!full tree level

      call FF_TX_(f3_type,f3_gen,f4_type,f4_gen,FF_TX)
      call FF_UX_(f3_type,f3_gen,f4_type,f4_gen,FF_UX)

      if(channel.eq.99) then

        AmpTree = FF_TX(6)+ FF_UX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = 0d0
          AmpTreeCh(i+5) = 0d0
          AmpTreeCh(i+10) = 0d0
          AmpTreeCh(i+15) = FF_TX(i)
          AmpTreeCh(i+20) = FF_UX(i)

        enddo

        AmpTree = AmpTreeCh(channel)

      endif

      AmpProp = 0d0
      AmpVert = 0d0
      AmpBox = 0d0
                                                                                                      
!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine FF_TX_(f3_type,f3_gen,f4_type,f4_gen,FF_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,f3_type,f3_gen,f4_type,f4_gen
      double precision colorfacST,Amp_TQ,Amp_TS,Amp_TV,Amp_TT,Amp_TU,FF_TX(6)

      Amp_TQ = 0d0
      Amp_TS = 0d0
      Amp_TV = 0d0
      Amp_TT = 0d0
      Amp_TU = 0d0

!==========================================================t*q=====================================================================

      Amp_TQ = 0d0    

!========================================================t*Scalar==================================================================

      Amp_TS = 0d0

!========================================================t*Vector==================================================================
 
      Amp_TV = 0d0

!========================================================t*t=======================================================================  
!Gluino exchange 

!tt(Gluino-exchange: fermion-fermion and antifermion-antifermion-final state)
      Amp_TT = Amp_TT + (2d0/3d0)*( dconjg(GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*
     &     (dconjg(GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &      (m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (MGl*(m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 
     &           2*m4*tman*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) - 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*
     &            GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + m4*MGl*(-m1**2 + m3**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))
     &      - MGl*dconjg(GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &      (m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))
     &         + GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))) + 
     &       dconjg(GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*
     &   (MGl*dconjg(GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &      (m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &     dconjg(GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &      (-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (2*m4*tman*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (m4*MGl*(m1**2 - m3**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 
     &           (-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*
     &            GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))))/((tman-MGl2)*(tman-MGl2))   


!Neutralino/Chargino-exchange

      do i= 1, 6
        do j=1, 6        
 
!tt(Neutralino/Chargino-exchange: fermion-fermion and antifermion-antifermion-final state)
          Amp_TT = Amp_TT + NC*(-(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &        (dconjg(AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &           (m3*(-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &             (-m1**2 - m2**2 + tman + uman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &          m4*dconjg(AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &           ((-m1**2 + m3**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &             2*m3*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j))) + 
     &       AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &        (dconjg(AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &           ((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*
     &              dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &             m3*(-m2**2 + m4**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &          m4*dconjg(AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &           (2*m3*tman*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &             (-m1**2 + m3**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j))))) - 
     &           AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (m4*dconjg(AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &         (2*m3*tman*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &           (-m1**2 + m3**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &         ((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*
     &            dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &           m3*(-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &      (m4*dconjg(AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &         ((-m1**2 + m3**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &           2*m3*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &         (m3*(-m2**2 + m4**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &           (-m1**2 - m2**2 + tman + uman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)))))/((tman-MNeuCha2(i))*(tman-MNeuCha2(j)))

        enddo
      enddo

!========================================================t*u=======================================================================
!Gluino exchange 

!tu+ut(Gluino-exchange: fermion-fermion and antifermion-antifermion-final state)
      Amp_TU = Amp_TU + (-2d0/9d0)*dreal(  (  dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (MGl*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &     dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (-(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (MGl*(-m1**2 + m4**2 + uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*MGl*(m2**2 - m3**2 - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))) + 
     &        dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (-(MGl*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &        (m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))
     &            + GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))) + 
     &     dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 + m4**2 + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (m4*MGl*(m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))))   )
     &        /((tman-MGl2)*(uman-MGl2))    )

!Neutralino/Chargino-exchange
!t/u(Neutralino/Chargino)*t/u(Gluino)=0 due to color-sums

      do i= 1, 6

!t(Neutralino/Chargino)*u(Gluino)  
        Amp_TU = Amp_TU + (4d0/3d0)*dreal(   (    AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (-(m4*dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &           (m3*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &             MGl*(-m1**2 + m3**2 + tman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &        dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (m3*MGl*(m2**2 - m4**2 - tman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         ((m2**2 - m3**2 - uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) - 
     &           2*m3*MGl*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &        dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           m3*(m1**2 - m4**2 - uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))))*MNeuCha(i)) + 
     &       AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (-(m4*dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &           (MGl*(-m1**2 + m3**2 + tman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &             m3*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &        dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           m3*MGl*(m2**2 - m4**2 - tman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) - 
     &     AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (2*m3*MGl*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           (-m2**2 + m3**2 + uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &        dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (m3*(-m1**2 + m4**2 + uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))))*MNeuCha(i))   ) 
     &        /((tman-MNeuCha2(i))*(uman-MGl2))    )

!t(Gluino)*u(Neutralino/Chargino) 
        Amp_TU = Amp_TU + (4d0/3d0)*dreal(  (  dconjg(AR(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*(dconjg(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &      (-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 + m4**2 + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*MGl*(m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &     dconjg(AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))*MNeuCha(i)) + 
     &       dconjg(AL(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*(dconjg(AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &      (-(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(-m1**2 + m4**2 + uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*MGl*(m2**2 - m3**2 - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) - 
     &     dconjg(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))*MNeuCha(i))  ) /((uman-MNeuCha2(i))*(tman-MGl2))  )

        do j=1, 6        
 
!tu(Neutralino/Chargino-exchange: fermion-fermion and antifermion-antifermion-final state)
          Amp_TU = Amp_TU + dreal( (   AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &      (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(m2**2 - m3**2 - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (m1**2 + m2**2 - tman - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2 - m4**2 - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) - 
     &           2*m4*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(-(m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &             (-m2**2 + m4**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           m4*(m1**2 - m3**2 - tman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)))) + 
     &  AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           m4*(m1**2 - m3**2 - tman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) - 
     &        m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &      (m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2 - m4**2 - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) - 
     &           2*m4*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(m2**2 - m3**2 - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (m1**2 + m2**2 - tman - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))))  )/((tman-MNeuCha2(i))*(uman-MNeuCha2(j)))  )

        enddo
      enddo

!==================================================================================================================================

      FF_TX(1) = Amp_TQ
      FF_TX(2) = Amp_TS
      FF_TX(3) = Amp_TV
      FF_TX(4) = Amp_TT
      FF_TX(5) = Amp_TU

      FF_TX(6) = FF_TX(1) + FF_TX(2) + FF_TX(3) + FF_TX(4) + FF_TX(5)
!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine FF_UX_(f3_type,f3_gen,f4_type,f4_gen,FF_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,f3_type,f3_gen,f4_type,f4_gen
      double precision colorfacST,Amp_UQ,Amp_US,Amp_UV,Amp_UT,Amp_UU,FF_UX(6)

      Amp_UQ = 0d0
      Amp_US = 0d0
      Amp_UV = 0d0
      Amp_UT = 0d0
      Amp_UU = 0d0

!==========================================================u*q=====================================================================

      Amp_UQ = 0d0    

!========================================================u*Scalar==================================================================

      Amp_US = 0d0

!========================================================u*Vector==================================================================
 
      Amp_UV = 0d0

!========================================================u*t=======================================================================  
!Gluino exchange 

!tu+ut(Gluino-exchange: fermion-fermion and antifermion-antifermion-final state)
      Amp_UT = Amp_UT + (-2d0/9d0)*dreal(  (  dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (MGl*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &     dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (-(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (MGl*(-m1**2 + m4**2 + uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*MGl*(m2**2 - m3**2 - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))) + 
     &     dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (-(MGl*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &        (m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))
     &            + GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))) + 
     &     dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &           (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 + m4**2 + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (m4*MGl*(m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))))   )
     &        /((tman-MGl2)*(uman-MGl2))    )

!Neutralino/Chargino-exchange
!t/u(Neutralino/Chargino)*t/u(Gluino)=0 due to color-sums

      do i= 1, 6

!t(Neutralino/Chargino)*u(Gluino)  
        Amp_UT = Amp_UT + (4d0/3d0)*dreal(   (    AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (-(m4*dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &           (m3*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &             MGl*(-m1**2 + m3**2 + tman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &        dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (m3*MGl*(m2**2 - m4**2 - tman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         ((m2**2 - m3**2 - uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) - 
     &           2*m3*MGl*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &        dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           m3*(m1**2 - m4**2 - uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))))*MNeuCha(i)) + 
     &     AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (-(m4*dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &           (MGl*(-m1**2 + m3**2 + tman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &             m3*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &        dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           m3*MGl*(m2**2 - m4**2 - tman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) - 
     &     AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (2*m3*MGl*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           (-m2**2 + m3**2 + uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &        dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &         (m3*(-m1**2 + m4**2 + uman)*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))))*MNeuCha(i))   ) 
     &        /((tman-MNeuCha2(i))*(uman-MGl2))    )

!t(Gluino)*u(Neutralino/Chargino) 
        Amp_UT = Amp_UT + (4d0/3d0)*dreal(  (  dconjg(AR(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*(dconjg(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &      (-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             MGl*(-m1**2 + m4**2 + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*MGl*(m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           (m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &     dconjg(AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         ((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))*MNeuCha(i)) + 
     &     dconjg(AL(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*(dconjg(AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &      (-(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(-m1**2 + m4**2 + uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &             m4*(-m1**2 - m2**2 + m3**2 + m4**2)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           m4*MGl*(m2**2 - m3**2 - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) - 
     &     dconjg(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &         (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))*MNeuCha(i))  ) /((uman-MNeuCha2(i))*(tman-MGl2))  )

        do j=1, 6        
 
!ut(Neutralino/Chargino-exchange: fermion-fermion and antifermion-antifermion-final state)
          Amp_UT = Amp_UT + dreal( (   AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &      (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(m2**2 - m3**2 - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (m1**2 + m2**2 - tman - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2 - m4**2 - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) - 
     &           2*m4*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(-(m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &             (-m2**2 + m4**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           m4*(m1**2 - m3**2 - tman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)))) + 
     &  AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &      (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2*m2**2 - m3**2*m4**2 - tman*uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           m4*(m1**2 - m3**2 - tman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) - 
     &        m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(-m1**2 - m2**2 + m3**2 + m4**2)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*MNeuCha(i)*
     &      (m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((m1**2 - m4**2 - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) - 
     &           2*m4*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(m2**2 - m3**2 - uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (m1**2 + m2**2 - tman - uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))))  )/((tman-MNeuCha2(i))*(uman-MNeuCha2(j)))  )

        enddo
      enddo

!========================================================u*u=======================================================================
!Gluino exchange 

      Amp_UU = Amp_UU + (2d0/3d0)*(     dconjg(GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (m4*GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         (MGl*(m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) - 
     &           2*m3*uman*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) - 
     &        GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         ((m1**2*(m2**2 - m3**2) - m2**2*m4**2 + m3**2*m4**2 + m3**2*uman + m4**2*uman - tman*uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) + 
     &           m3*MGl*(-m1**2 + m4**2 + uman)*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) - 
     &     MGl*dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (m4*GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         (2*m3*MGl*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) + (-m2**2 + m3**2 + uman)*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))
     &         + GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         (m3*(-m1**2 + m4**2 + uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) + 
     &           MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)))) + 
     &    dconjg(GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &   (MGl*dconjg(GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (m4*GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         ((m2**2 - m3**2 - uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) - 2*m3*MGl*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)) + 
     &        GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         (MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) + 
     &           m3*(m1**2 - m4**2 - uman)*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &     dconjg(GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &      (-(m4*GluR(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &           (2*m3*uman*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) + 
     &             MGl*(-m2**2 + m3**2 + uman)*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))) + 
     &        GluL(sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &         (m3*MGl*(m1**2 - m4**2 - uman)*GluL(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen) - 
     &           (m1**2*(m2**2 - m3**2) - m2**2*m4**2 + m3**2*m4**2 + m3**2*uman + m4**2*uman - tman*uman )*GluR(sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))))   )
     &        /((uman-MGl2)*(uman-MGl2))    

!Neutralino/Chargino-exchange

      do i= 1, 6
        do j=1, 6        
 
!uu(Neutralino/Chargino-exchange: fermion-fermion and antifermion-antifermion-final state)
          Amp_UU = Amp_UU + NC*( -(AL(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)*(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*MNeuCha(i)*
     &        (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           (m4*(-m2**2 + m3**2 + uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &             (-m1**2 - m2**2 + tman + uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &          m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           ((-m1**2 + m4**2 + uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &             2*m4*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &       AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &        (dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           (dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &              (m1**2*(m2**2 - m3**2) - m2**2*m4**2 + m3**2*m4**2 + m3**2*uman + m4**2*uman - tman*uman ) + 
     &             m4*(-m2**2 + m3**2 + uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &          m3*dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &           (2*m4*uman*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &             (-m1**2 + m4**2 + uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))))) - 
     &         AR(i,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen)*(AL(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*
     &      (m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (2*m4*uman*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (-m1**2 + m4**2 + uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*
     &            (m1**2*(m2**2 - m3**2) - m2**2*m4**2 + m3**2*m4**2 + m3**2*uman + m4**2*uman - tman*uman ) + 
     &           m4*(-m2**2 + m3**2 + uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))) + 
     &     AR(i,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)*MNeuCha(i)*
     &      (m3*dconjg(AL(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         ((-m1**2 + m4**2 + uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           2*m4*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j)) + 
     &        dconjg(AR(j,sfmeigen2,sftype2,f3_type,sfgen2,f3_gen))*
     &         (m4*(-m2**2 + m3**2 + uman)*dconjg(AL(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen)) + 
     &           (-m1**2 - m2**2 + tman + uman)*dconjg(AR(j,sfmeigen1,sftype1,f4_type,sfgen1,f4_gen))*MNeuCha(j))))   )/((uman-MNeuCha2(i))*(uman-MNeuCha2(j))) 

        enddo
      enddo

!==================================================================================================================================

      FF_UX(1) = Amp_UQ
      FF_UX(2) = Amp_US
      FF_UX(3) = Amp_UV
      FF_UX(4) = Amp_UT
      FF_UX(5) = Amp_UU

      FF_UX(6) = FF_UX(1) + FF_UX(2) + FF_UX(3) + FF_UX(4) + FF_UX(5)

!==================================================================================================================================

      end