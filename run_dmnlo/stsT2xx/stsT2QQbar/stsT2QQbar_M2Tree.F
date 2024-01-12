      double precision function stsT2QQbar_M2Tree()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double complex DenSV
      double precision colfac

! gluon propagator
      DenSV = sman


! --------------------- Amplitude squared for S-channel with gluon -----------

      colfac=(Nc**2-1)/4d0

      stsT2QQbar_M2Tree= - colfac*((sqrtS**4*((-GS)*      ((4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred  + 4*tred**2)*
     &            (-GS) ) + (-GS)*((4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred + 4*tred**2)*
     &           (-GS)))*dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     &      SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/ (dconjg(DenSV)*DenSV))

	end

      double precision function stsT2QQbar_M2TreeDipole(sdip,tdip,udip)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
      double precision sdip,tdip,udip

      stsT2QQbar_M2TreeDipole = -((Nc*CF*(4*m1**2*sdip - sdip**2 + (tdip - udip)**2)*gs**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/sdip**2)

      end

! total cross section --  octet (singlet is zero)
      double precision function stsT2QQbar_TreeTotalOctet()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
      
      double precision rho,beta

      rho = (4*m1**2)/sman
      beta = dsqrt(1-rho)

      stsT2QQbar_TreeTotalOctet = (8*Pi*Alfas2*beta)/(27*sman)

      end