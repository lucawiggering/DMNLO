	function stst2gg_M2Tree()	
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision M2,Krondelta
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      double complex stst2gg_M2Tree
      integer i,j,k  
 
      M2=MSf2(iSf1,itt1,igen1)
      !Propagators

      DenS = sman 

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo

      
!!This document contains the triple checked tree-level, 
!!based on calculations by SS, JB and LW. In RenScheme99,
!!the micrOMEGAs result is exactly the same.

! ! !---------------------|Ms|^2-----------
! 
      stst2gg_M2Tree = (CF*CV*GS**4*Nc*(2*(-9 + 5*FPSk1)*M2**2 - 8*(-2 + FPSk1)*M2*sman - 8*sman**2 + 4*FPSk1*sman**2 - 2*M2*((-9 + 5*FPSk1)*M2 + 4*(-2 + FPSk1)*sman) + 5*tman**2 + 
     -      3*FPSk1*tman**2 - 10*tman*uman - 6*FPSk1*tman*uman +5*uman**2 + 3*FPSk1*uman**2)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (2d0*DenS*dconjg(DenS))
! ! ---------------------|Mt|^2-----------

      stst2gg_M2Tree = stst2gg_M2Tree+(GS**4*(-1 + Nc**2 +CF*Nc*(-4 + CV*Nc + Nc**2))*(M2 + tman)*(sman*(M2 + tman) - FPSk1*(M2 - tman)*(2*M2 - sman - 2*uman))*Krondelta(igen1,igen2)*
     -    Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(Nc*sman*dconjg(DenSft(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1))

! ! ! -------------------------------------|Mu|^2------------------------------------------------

      stst2gg_M2Tree = stst2gg_M2Tree+(GS**4*(-1 + Nc**2 +CF*Nc*(-4 + CV*Nc + Nc**2))*(M2 + uman)*(-(FPSk1*(2*M2 - sman - 2*tman)*(M2 - uman)) + sman*(M2 + uman))*Krondelta(igen1,igen2)*
     -    Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(Nc*sman*dconjg(DenSfu(isf1,itt1,igen1))*DenSfu(isf1,itt1,igen1))
! 
! !---------------------|Mv|^2-----------

      stst2gg_M2Tree = stst2gg_M2Tree+(-2*(-2 + FPSk1)*(-1 + Nc**2 + CF*Nc*(-4 + Nc**2))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/Nc
     

!--------------------- 2 Re(Ms Mt^\ast)-----------

      stst2gg_M2Tree = stst2gg_M2Tree-Real((CF*CV*GS**4*Nc*(sman*(-18*M2**2 - 4*sman**2 - 7*sman*tman + 2*tman**2 - 5*sman*uman - 2*tman*uman + M2*(14*sman + 11*tman + 13*uman) + 
     -          M2*(-30*M2 + 14*sman + 11*tman + 13*uman)) + FPSk1*(2*sqrtS**6 + sman**2*tman + 6*tman**3 + M2**2*(12*M2 + 5*sman - 4*tman - 8*uman) + 
     -          5*sman**2*uman - 2*sman*tman*uman - 4*tman**2*uman + 2*sman*uman**2 - 2*tman*uman**2 - 
     -          M2*(6*sman**2 + 7*sman*tman + 12*tman**2 + 15*sman*uman - 6*tman*uman + 6*uman**2) + 
     -          M2*(-12*M2**2 -8*sman**2 + sman*tman - 2*tman**2 + 2*M2*(8*sman + 3*tman - 3*uman) - 3*sman*uman + 6*tman*uman + 8*uman**2) + 
     -          3*M2**2*(sman + 2*(tman + uman))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(2d0*DenS*sman*dconjg(DenSft(isf1,itt1,igen1))))

! ! ! !---------------------2 Re(Ms Mu^\ast)-----------

      stst2gg_M2Tree = stst2gg_M2Tree-Real((CF*CV*GS**4*Nc*(sman*(-18*M2**2 - 4*sman**2 - 5*sman*tman - 7*sman*uman - 2*tman*uman + 2*uman**2 + M2*(14*sman + 13*tman + 11*uman) + 
     -          M2*(-30*M2 + 14*sman + 13*tman + 11*uman)) + FPSk1*(2*sqrtS**6 + 5*sman**2*tman + 2*sman*tman**2 + sman**2*uman - 2*sman*tman*uman - 
     -          2*tman**2*uman - 4*tman*uman**2 + 6*uman**3 + M2*(-8*sman**2 + 8*tman**2 + 6*tman*uman - 2*uman**2 + sman*(-3*tman + uman)) + 
     -          3*M2**2*(-4*M2 + sman + 2*(tman + uman)) + M2**2*(5*sman-4*(2*tman + uman)) + 
     -          M2*(12*M2**2 - 6*sman**2 - 15*sman*tman -6*tman**2 - 7*sman*uman + 6*tman*uman - 12*uman**2 + 2*M2*(8*sman - 3*tman + 3*uman))))*
     -     Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(2d0*DenS*sman*dconjg(DenSfu(isf1,itt1,igen1))))
!  
! ! !---------------------2 Re(Mt Mu^\ast)-----------
 	
      stst2gg_M2Tree = stst2gg_M2Tree+Real((GS**4*(-(1/Nc) + Nc + CF*(-4 - CV*Nc + Nc**2))*(-6*M2 + 2*sman + tman + uman)*
     -    (sman*(-6*M2 + 2*sman + tman + uman) + FPSk1*(-(sman*tman) - 2*tman**2 - sman*uman - 2*uman**2 + M2*(sman + 2*(tman + uman)) + 
     -         M2*(-4*M2 +sman + 2*(tman + uman))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (2d0*sman*dconjg(DenSfu(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)))

! ! ! --------------------- 2 Re(Mt Mv^\ast)-----------

      stst2gg_M2Tree = stst2gg_M2Tree -Real((GS**2*(-1 + Nc**2 + CF*Nc*(-4 + Nc**2))*(sman*(12*M2 - 5*sman - 4*uman) + 
     -       FPSk1*(sman**2 + 4*tman**2 + 4*sman*uman + 4*uman**2 - 2*M2*(sman + 2*(tman + uman)) + M2*(8*M2 - 2*(sman + 2*(tman + uman)))))*
     -     Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2d0*Nc*sman*DenSft(isf1,itt1,igen1)))
      
! ! ---------------------2 Re(Mu Mv^\ast)-----------

      stst2gg_M2Tree = stst2gg_M2Tree-Real((GS**2*(-1 + Nc**2 + CF*Nc*(-4 +Nc**2))*(sman*(12*M2 - 5*sman - 4*tman) + 
     -       FPSk1*(sman**2 + 4*sman*tman + 4*tman**2 + 4*uman**2 - 2*M2*(sman + 2*(tman + uman)) + M2*(8*M2 - 2*(sman + 2*(tman + uman)))))*
     -     Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2d0*Nc*sman*DenSfu(isf1,itt1,igen1)))

!! --------------------   -|Mgh1|^2    -------------------------------
      stst2gg_M2Tree = stst2gg_M2Tree -(1d0-FPSk1)*(CF*CV*GS**4*Nc*(-tman + uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS))

!! --------------------   -|Mgh2|^2    -------------------------------
      stst2gg_M2Tree = stst2gg_M2Tree -(1d0-FPSk1)*(CF*CV*GS**4*Nc*(tman- uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS))
      
          
      end


! ! 	Color decomposition

      double precision function stst2gg_M2TreeSinglet()	
      implicit none  
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision M2,Krondelta
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      integer i,j,k  

      M2=MSf2(iSf1,itt1,igen1)
      !Propagators

      DenS = sman 

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo

!---------------- |(Mt+Mu+Mv)^{(1)}|^2   -------------------------- 

      stst2gg_M2TreeSinglet=  dble(((-1 + Nc**2)*(GS**2*dconjg(DenSfu(isf1,itt1,igen1))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*
     -       (4*GS**2*(M2 + tman)**2*DenSfu(isf1,itt1,igen1) + DenSft(isf1,itt1,igen1)*
     -          (GS**2*(-6*M2 + 2*sman + tman + uman)**2 +(-12*M2 + 5*sman + 4*uman)*DenSfu(isf1,itt1,igen1)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))) + 
     -      dconjg(DenSft(isf1,itt1,igen1))*(GS**2*(GS**2*(-6*M2 +2*sman+ tman + uman)**2 + 
     -            (-12*M2 + 5*sman +4*uman)*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)))*DenSfu(isf1,itt1,igen1)*
     -          Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2) + 
     -         DenSft(isf1,itt1,igen1)*(16*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))*DenSfu(isf1,itt1,igen1)*
     -             SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2) + GS**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*
     -             (4*GS**2*(M2 + uman)**2 + (-12*M2 + 5*sman + 4*tman)*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)) + 
     -               (-12*M2 + 5*sman + 4*tman)*DenSfu(isf1,itt1,igen1)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))))))/
     -  (4d0*Nc*dconjg(DenSft(isf1,itt1,igen1))*dconjg(DenSfu(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)*DenSfu(isf1,itt1,igen1)))

      end

      double precision function stst2gg_M2TreeSOctet()	
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
 
      double precision M2,Krondelta
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      integer i,j,k  

      M2=MSf2(iSf1,itt1,igen1)
      !Propagators

      DenS = sman 

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo


!---------------- |(Mt+Mu+Mv)^{(8_S)}|^2   --------------------------  

      stst2gg_M2TreeSOctet= dble((CF*(-4 + Nc**2)*(GS**2*dconjg(DenSfu(isf1,itt1,igen1))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*
     -       (4*GS**2*(M2 + tman)**2*DenSfu(isf1,itt1,igen1) + DenSft(isf1,itt1,igen1)*
     -          (GS**2*(-6*M2 + 2*sman + tman + uman)**2 + (-12*M2+ 5*sman + 4*uman)*DenSfu(isf1,itt1,igen1)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))) + 
     -      dconjg(DenSft(isf1,itt1,igen1))*(GS**2*(GS**2*(-6*M2 + 2*sman + tman + uman)**2 + 
     -            (-12*M2 + 5*sman +4*uman)*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)))*DenSfu(isf1,itt1,igen1)*
     -          Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2) + 
     -         DenSft(isf1,itt1,igen1)*(16*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))*DenSfu(isf1,itt1,igen1)*
     -             SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2) + GS**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2)*
     -             (4*GS**2*(M2 + uman)**2 + (-12*M2 + 5*sman + 4*tman)*dconjg(DenSfu(isf1,itt1,igen1))*dconjg(SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)) + 
     -               (-12*M2 + 5*sman + 4*tman)*DenSfu(isf1,itt1,igen1)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))))))/
     -  (4d0*dconjg(DenSft(isf1,itt1,igen1))*dconjg(DenSfu(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)*DenSfu(isf1,itt1,igen1)))

      end

      double precision function stst2gg_M2TreeAOctet()	
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision M2,Krondelta
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      integer i,j,k  

      M2=MSf2(iSf1,itt1,igen1)
      !Propagators

      DenS = sman 

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo

!---------------- |(Ms+Mt+Mu)^{(8_A)}|^2   --------------------------      

      stst2gg_M2TreeAOctet = dble(-(CF*CV*GS**4*Nc*(dconjg(DenSft(isf1,itt1,igen1))*dconjg(DenSfu(isf1,itt1,igen1))*
     -        (DenS*(sman*(-18*M2**2 - 4*sman**2 - 7*sman*tman + 2*tman**2 - 5*sman*uman - 2*tman*uman + M2*(14*sman + 11*tman + 13*uman) + 
     -                M2*(-30*M2 + 14*sman + 11*tman + 13*uman)) + 
     -             FPSk1*(2*sqrtS**6 + sman**2*tman + 6*tman**3 + M2**2*(12*M2 + 5*sman - 4*tman - 8*uman) + 5*sman**2*uman - 2*sman*tman*uman - 4*tman**2*uman + 
     -                2*sman*uman**2 - 2*tman*uman**2 - M2*(6*sman**2 + 7*sman*tman + 12*tman**2 + 15*sman*uman - 6*tman*uman + 6*uman**2) + 
     -                M2*(-12*M2**2 - 8*sman**2 + sman*tman - 2*tman**2 + 2*M2*(8*sman + 3*tman - 3*uman) - 3*sman*uman + 6*tman*uman + 8*uman**2) + 
     -                3*M2**2*(sman+ 2*(tman+ uman))))*DenSfu(isf1,itt1,igen1) + 
     -          DenSft(isf1,itt1,igen1)*(DenS*(sman*(-18*M2**2 - 4*sman**2 - 5*sman*tman - 7*sman*uman - 2*tman*uman + 2*uman**2 + 
     -                   M2*(14*sman + 13*tman + 11*uman) + M2*(-30*M2+ 14*sman + 13*tman + 11*uman)) + 
     -                FPSk1*(2*sqrtS**6 + 5*sman**2*tman + 2*sman*tman**2 + sman**2*uman - 2*sman*tman*uman - 2*tman**2*uman - 4*tman*uman**2 + 6*uman**3 + 
     -                   M2*(-8*sman**2 + 8*tman**2 + 6*tman*uman - 2*uman**2 + sman*(-3*tman + uman)) + 3*M2**2*(-4*M2 + sman + 2*(tman + uman)) + 
     -                   M2**2*(5*sman - 4*(2*tman + uman)) + M2*(12*M2**2 - 6*sman**2 - 15*sman*tman - 6*tman**2 - 7*sman*uman + 6*tman*uman - 12*uman**2 + 
     -                      2*M2*(8*sman - 3*tman + 3*uman)))) - 2*sman*
     -              (2*(-9 + 5*FPSk1)*M2**2 - 8*(-2 + FPSk1)*M2*sman - 8*sman**2 + 4*FPSk1*sman**2 - 2*M2*((-9 + 5*FPSk1)*M2 + 4*(-2 + FPSk1)*sman) + 5*tman**2 + 
     -                3*FPSk1*tman**2 - 10*tman*uman - 6*FPSk1*tman*uman + 5*uman**2 + 3*FPSk1*uman**2)*DenSfu(isf1,itt1,igen1))) + 
     -       dconjg(DenS)*(dconjg(DenSfu(isf1,itt1,igen1))*(-4*DenS*(M2 + tman)*(sman*(M2 + tman) - FPSk1*(M2 - tman)*(2*M2 - sman - 2*uman))*
     -              DenSfu(isf1,itt1,igen1) + DenSft(isf1,itt1,igen1)*
     -              (DenS*(6*M2 - 2*sman - tman - uman)*(-(sman*(-6*M2 +2*sman + tman + uman)) + 
     -                   FPSk1*(sman*tman + 2*tman**2 + sman*uman + 2*uman**2 + M2*(4*M2 - sman - 2*(tman + uman)) - M2*(sman + 2*(tman + uman)))) + 
     -                (sman*(-18*M2**2 - 4*sman**2 - 7*sman*tman + 2*tman**2 - 5*sman*uman - 2*tman*uman + M2*(14*sman + 11*tman + 13*uman) + 
     -                      M2*(-30*M2 + 14*sman + 11*tman + 13*uman)) + 
     -                   FPSk1*(2*sqrtS**6 + sman**2*tman + 6*tman**3 + M2**2*(12*M2 + 5*sman - 4*tman - 8*uman) + 5*sman**2*uman - 2*sman*tman*uman - 
     -                      4*tman**2*uman + 2*sman*uman**2 - 2*tman*uman**2 - 
     -                      M2*(6*sman**2 + 7*sman*tman + 12*tman**2 + 15*sman*uman - 6*tman*uman + 6*uman**2) + 
     -                      M2*(-12*M2**2 - 8*sman**2 + sman*tman- 2*tman**2 + 2*M2*(8*sman + 3*tman - 3*uman) - 3*sman*uman + 6*tman*uman + 8*uman**2) + 
     -                      3*M2**2*(sman + 2*(tman + uman))))*DenSfu(isf1,itt1,igen1))) + 
     -          dconjg(DenSft(isf1,itt1,igen1))*(DenS*(6*M2 - 2*sman- tman - uman)*
     -              (-(sman*(-6*M2 + 2*sman + tman + uman)) + 
     -                FPSk1*(sman*tman + 2*tman**2 + sman*uman + 2*uman**2+ M2*(4*M2 - sman - 2*(tman + uman)) - M2*(sman + 2*(tman + uman))))*
     -              DenSfu(isf1,itt1,igen1) + DenSft(isf1,itt1,igen1)*
     -              (-4*DenS*(M2 + uman)*(-(FPSk1*(2*M2 - sman - 2*tman)*(M2 - uman)) + sman*(M2 + uman)) + 
     -                (sman*(-18*M2**2 - 4*sman**2 - 5*sman*tman - 7*sman*uman - 2*tman*uman + 2*uman**2 + M2*(14*sman + 13*tman + 11*uman) + 
     -                      M2*(-30*M2 + 14*sman + 13*tman + 11*uman)) + 
     -                   FPSk1*(2*sqrtS**6 + 5*sman**2*tman + 2*sman*tman**2 + sman**2*uman - 2*sman*tman*uman - 2*tman**2*uman - 4*tman*uman**2 + 6*uman**3 + 
     -                      M2*(-8*sman**2 + 8*tman**2 + 6*tman*uman - 2*uman**2 + sman*(-3*tman + uman)) + 3*M2**2*(-4*M2 + sman + 2*(tman + uman)) + 
     -                      M2**2*(5*sman - 4*(2*tman + uman)) + M2*
     -                       (12*M2**2 - 6*sman**2 - 15*sman*tman- 6*tman**2 - 7*sman*uman + 6*tman*uman - 12*uman**2 + 2*M2*(8*sman - 3*tman + 3*uman))))*
     -                 DenSfu(isf1,itt1,igen1)))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4d0*DenS*sman*dconjg(DenS)*dconjg(DenSft(isf1,itt1,igen1))*dconjg(DenSfu(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)*DenSfu(isf1,itt1,igen1)))

!! --------------------   -|Mgh1|^2    -------------------------------
      stst2gg_M2TreeAOctet = stst2gg_M2TreeAOctet -(1d0-FPSk1)*dble((CF*CV*GS**4*Nc*(-tman + uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS)))

!! --------------------   -|Mgh2|^2    -------------------------------
      stst2gg_M2TreeAOctet = stst2gg_M2TreeAOctet -(1d0-FPSk1)*dble((CF*CV*GS**4*Nc*(tman- uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS)))
              
      end

! total cross section -- singlet
      double precision function stst2gg_TreeTotalSinglet()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision rho,beta,M2

      M2=MSf2(iSf1,itt1,igen1)

      rho = 4*MSf2(iSf1,itt1,igen1)/sman
      beta = dsqrt(1-rho)

      stst2gg_TreeTotalSinglet = 8*Pi*Alfas2/(27*sman*beta**2)*(beta*(1+rho)+rho*(rho-2)*atanh(beta))

      end

! total cross section -- symmetric octet
      double precision function stst2gg_TreeTotalSOctet()
      implicit none

      double precision stst2gg_TreeTotalSinglet

      stst2gg_TreeTotalSOctet = 5*stst2gg_TreeTotalSinglet()/2d0
      
      end

! total cross section -- antisymmetric octet
      double precision function stst2gg_TreeTotalAOctet()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
      
      double precision rho,beta

      rho = 4*MSf2(iSf1,itt1,igen1)/sman
      beta = dsqrt(1-rho)

      stst2gg_TreeTotalAOctet = 4*Pi*Alfas2/(9*sman*beta**2)*(beta*(1+8*rho)-3*rho*(rho+2)*atanh(beta))

      end


      double precision function stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
      double precision sdip,tdip,udip,M2

      M2=MSf2(iSf1,itt1,igen1)

 !************************************************   Ms^ast  Ms   ************************************************

      stsT2gg_M2TreeDipoleT1T3 = ((GGG**2*Nc**2*(-1 + Nc**2)*(-4*M2*sdip 
     &   + sdip**2 - 2*(tdip - udip)**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*sdip**2))

!************************************************   2 Re(Ms Mt^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble(-0.125*(GGG*Nc**2*(-1 + Nc**2)*(sdip
     &   **3 + 3*sdip**2*tdip - sdip*(tdip - udip)**2 + 4*M2**2*(3*sdip - tdip + udip) + tdip*(-3*tdip**2 
     &   + 2*tdip*udip + udip**2) - M2*(7*sdip**2 - 7*tdip**2 + 6*tdip*udip + udip**2 + 4*sdip*(2*tdip + udip)))
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble(0)

!************************************************   2 Re(Ms Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble((GGG*Nc**2*(-1 + Nc**2)*(tdip - udip)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*sdip))

!************************************************   Mt^ast  Mt   ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 +(((-1 + Nc**2)**2*(M2 + tdip)*(M2**2 + tdip
     &   *udip - M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble(((-1 + Nc**2)*(-6*M2 + 2*sdip + tdip 
     &   + udip)*(-2*M2**2 + sdip**2 - tdip**2 - udip**2 + 2*M2*(-sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(4*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble(-0.25d0*((2 - 3*Nc**2 + Nc**4)*(2*M2
     &   **2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(sdip*(M2 - tdip)))

!************************************************   Mu^ast  Mu   ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 +(((-1 + Nc**2)*(M2 + udip)*(-M2**2 - tdip
     &   *udip + M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 + 2*dble(((-1 + Nc**2)*(2*M2**2 - sdip**2 
     &   + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2
     &   *sdip*(M2 - udip)))

!************************************************   Mv^ast  Mv   ************************************************

      stsT2gg_M2TreeDipoleT1T3 = stsT2gg_M2TreeDipoleT1T3 +(-0.25d0*((4 - 5*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2))

      end



      double precision function stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2 

      M2=MSf2(iSf1,itt1,igen1)

!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2gg_M2TreeDipoleT2T3 = (-0.25d0*Nc*(GGG**2*Nc*(-1 
     &   + Nc**2)*(tdip - udip)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/sdip**2)

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble(-0.25d0*(GGG*Nc**2*(
     &   -1 + Nc**2)*(tdip - udip)*(2*M2**2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble((GGG*Nc**2*(-1 + Nc
     &   **2)*(-tdip + udip)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*sdip))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 +((2*(-1 + Nc**2)*(M2**2 
     &   + tdip*udip - M2*(sdip + tdip + udip))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(sdip**2*(M2  - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble(((-1 + Nc**2)*(2
     &   *M2**2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(2*sdip**2*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble(((-1 + Nc**2)*(2
     &   *M2**2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 +((-2*(-1 + Nc**2)**2*(M2
     &   **2 + tdip*udip - M2*(sdip + tdip + udip))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(sdip**2
     &   *(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 + 2*dble(-0.25d0*((2 - 3*Nc
     &   **2 + Nc**4)*(2*M2**2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2gg_M2TreeDipoleT2T3 = stsT2gg_M2TreeDipoleT2T3 +(-0.25d0*((4 - 5*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2))

      end

      double precision function stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2 

      M2=MSf2(iSf1,itt1,igen1)

!************************************************   Ms^ast  Ms   ************************************************

      stsT2gg_M2TreeDipoleT1T2 = ((GGG**2*(-1 + Nc**2)*(4*M2*sdip - sdip*
     &   *2 + 2*(tdip - udip)**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(2*sdip**2))

!************************************************   2 Re(Ms Mt^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble((GGG*(-1 + Nc**2)*(sdip**3 + 3*sdip
     &   **2*tdip - sdip*(tdip - udip)**2 + 4*M2**2*(3*sdip - tdip + udip) + tdip*(-3*tdip**2 + 2*tdip*udip 
     &   + udip**2) - M2*(7*sdip**2 - 7*tdip**2 + 6*tdip*udip + udip**2 + 4*sdip*(2*tdip + udip)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(8*sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble((GGG*(-1 + Nc**2)*(sdip**3 - sdip
     &   *(tdip - udip)**2 + 4*M2**2*(3*sdip + tdip - udip) + 3*sdip**2*udip + udip*(tdip**2 + 2*tdip*udip - 3
     &   *udip**2) - M2*(7*sdip**2 + tdip**2 + 6*tdip*udip - 7*udip**2 + 4*sdip*(tdip + 2*udip)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(8*sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble(0)

!************************************************   Mt^ast  Mt   ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 +(((-1 + Nc**2)*(M2 + tdip)*(M2**2 + tdip
     &   *udip - M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc**2*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Mu^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble(-0.25d0*((-1 + Nc**4)*(-6*M2 + 2*sdip 
     &   + tdip + udip)*(-2*M2**2 + sdip**2 - tdip**2 - udip**2 + 2*M2*(-sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(Nc**2*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble(-0.25d0*((-2 + Nc**2 + Nc**4)*(2*M2
     &   **2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(Nc**2*sdip*(M2 - tdip)))

!************************************************   Mu^ast  Mu   ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 +(((-1 + Nc**2)*(M2 + udip)*(M2**2 + tdip
     &   *udip - M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc**2*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Mv^ast)     ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 + 2*dble(-0.25d0*((-2 + Nc**2 + Nc**4)*(2*M2
     &   **2 - sdip**2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(Nc**2*sdip*(M2 - udip)))

!************************************************   Mv^ast  Mv   ************************************************

      stsT2gg_M2TreeDipoleT1T2 = stsT2gg_M2TreeDipoleT1T2 +(-0.5d0*((-2 + Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/Nc**2)

      end

      double precision function stsT2gg_M2TreeDipole(sdip,tdip,udip)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2 

      M2=MSf2(iSf1,itt1,igen1)

!************************************************   Ms^ast  Ms   ************************************************

      stsT2gg_M2TreeDipole = ((GGG**2*Nc*(-1 + Nc**2)*(4*M2*sdip - sdip**2 + 2
     &   *(tdip - udip)**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/sdip**2)

!************************************************   2 Re(Ms Mt^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble((GGG*Nc*(-1 + Nc**2)*(sdip**3 + 3*sdip**2
     &   *tdip - sdip*(tdip - udip)**2 + 4*M2**2*(3*sdip - tdip + udip) + tdip*(-3*tdip**2 + 2*tdip*udip + udip
     &   **2) - M2*(7*sdip**2 - 7*tdip**2 + 6*tdip*udip + udip**2 + 4*sdip*(2*tdip + udip)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(4*sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Mu^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble((GGG*Nc*(-1 + Nc**2)*(sdip**3 - sdip*(tdip 
     &   - udip)**2 + 4*M2**2*(3*sdip + tdip - udip) + 3*sdip**2*udip + udip*(tdip**2 + 2*tdip*udip - 3*udip
     &   **2) - M2*(7*sdip**2 + tdip**2 + 6*tdip*udip - 7*udip**2 + 4*sdip*(tdip + 2*udip)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(4*sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Mv^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble(0)

!************************************************   Mt^ast  Mt   ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole +((2*(-1 + Nc**2)**2*(M2 + tdip)*(-M2**2 - tdip*udip 
     &   + M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Mu^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble(-0.5d0*((-1 + Nc**2)*(-6*M2 + 2*sdip + tdip 
     &   + udip)*(-2*M2**2 + sdip**2 - tdip**2 - udip**2 + 2*M2*(-sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Mv^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble(((-2 + Nc**2)*(-1 + Nc**2)*(2*M2**2 - sdip
     &   **2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2
     &   *Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast  Mu   ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole +((2*(-1 + Nc**2)**2*(M2 + udip)*(-M2**2 - tdip*udip 
     &   + M2*(sdip + tdip + udip))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Mv^ast)     ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole + 2*dble(((-2 + Nc**2)*(-1 + Nc**2)*(2*M2**2 - sdip
     &   **2 + tdip**2 + 2*M2*(sdip - tdip - udip) + udip**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2.
     &   *Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast  Mv   ************************************************

      stsT2gg_M2TreeDipole = stsT2gg_M2TreeDipole +(((-2 + Nc**2)*(-1 + Nc**2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/Nc)

      end