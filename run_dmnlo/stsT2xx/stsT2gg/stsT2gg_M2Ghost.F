      subroutine stst2gg_calc_Aghost(A1,A2,k1k1,k2k2,k1k2)
      implicit none
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"
      double complex A1,A2
      double precision k1k1,k2k2,k1k2,k1Mk22,color

      ! k1Mk22 -> (k1-k2)^2=k1^2+k2^2-2*k1k2 
      k1Mk22=k1k1+k2k2-2d0*k1k2
      color=Nc/2d0

      call myCget(k2k2,k1Mk22,k1k1,0d0,0d0,0d0)
      call myBget(k1Mk22,0d0,0d0) 

! Gluon exchange 
      A1 = color*(-GS**3*(k1k2*(xc2 + xc22) + k2k2*xc12))/(16*Pi**2)
      A2 = color*(-GS**3*(xc00 + k1k2*(xc0 + xc1 + xc12 + xc2) + k2k2*(xc1 + xc11)))/(16*Pi**2)

! Ghost exchange
      A1 = A1 + color*(GS**3*((-xc2 + xc22)*k1k2 + (xc1 + xc12 + xc2)*k2k2))/(16*Pi**2)
      A2 = A2 + color*(GS**3*(-xb0 + xc00 + xc2*k1k1 - (xc1 - xc12 + xc2)*k1k2 + xc11*k2k2))/(16*Pi**2)

      end
      
      
      subroutine stst2gg_calc_M2VertexGh(M2Color)
	implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"

      double complex A1,A2,S1NLO,S2NLO,S1LO,S2LO
      double precision k1k1,k2k2,Krondelta,color,M2,M2Color(3)

      M2Color = 0d0
      M2=MSf2(iSf1,itt1,igen1)

! Scalar products for ghost-gluon vertex
      k1k1=0d0
      k2k2=0d0
      k1k2=-sman/2d0

! Calculation of partial amplitudes
      call stst2gg_calc_Aghost(A1,A2,k1k1,k2k2,k1k2)

      S1NLO = ((A1 + A2)*cI*GS*(tman - uman))/(2*sman)
      S2NLO = ((A1 + A2)*cI*GS*(tman - uman))/(2*sman)


! Leading order ghost amplitudes
      S1LO=(cI*GS**2*(uman -tman)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(2*sman)
      S2LO=(cI*GS**2*(uman -tman)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(2*sman)


! Interferences between tree level and NLO
      color=Nc*CV*CF   ! given by the scalar product <c^{8A}|c^{8A}>

! Interference of S1LO with S1NLO
      M2Color(3) = 2*dble(color*dconjg(S1LO)*S1NLO)
 
! Interference of S2LO with S2NLO
      M2Color(3) = M2Color(3) + 2*dble(color*dconjg(S2LO)*S2NLO)
      
      end

      subroutine stst2gg_calc_M2GhostBoxes(M2Color)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

! Since the ghost tree level amplitudes live in the 8A color subspace, only these components of the NLO amplitudes are considered.
      double complex S1NLO,S2NLO,S1LO,S2LO
      double precision color,M2,Krondelta,M2Color(3)

      M2Color = 0d0
      M2=MSf2(iSf1,itt1,igen1)

! ---------------------  Ghost exchanges right ------------------------
! Do not contribute (S1LO,S2LO only contribute to 8A whereas the "ghost exchanges" right live in singlet and symmetric octet)

! ---------------------  Ghost boxes -----------------------------------
      call LTDget(M2,0d0,0d0,M2,tman,sman,M2,0d0,0d0,0d0)
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
      color=-Nc/4d0 ! only 8A

! Ghost box 1
      S1NLO = cI*color*(-(GS**4*(-2*sman*xc1 - 4*M2*xc2 + M2**2*xd0 - 2*sman*xd00 - M2**2*xd11 - 2*M2**2*xd12 + M2*sman*xd12 + 2*M2**2*xd13 - M2*sman*xd2 - M2**2*xd22 + 
     &       M2*sman*xd22 + 2*M2**2*xd23 - M2*sman*xd23 + uman**2*(xd13 + xd23) + tman**2*(xd0 - xd1 + xd13 - xd2 + xd23 - xd3) - M2**2*xd33 + 
     &       uman*(2*xc0 + 2*xc2 - M2*xd1 + M2*xd11 + 2*M2*xd12 - 2*M2*xd13 - M2*xd2 + M2*xd22 - 2*M2*xd23 + sman*xd23 - M2*xd3 + M2*xd33) + 
     &       tman*(-2*xc0 + 2*xc2 - 2*M2*xd0 + M2*xd1 + uman*xd1 + M2*xd11 - uman*xd11 + 2*M2*xd12 - sman*xd12 - 2*uman*xd12 - 2*M2*xd13 + M2*xd2 + sman*xd2 + 
     &          uman*xd2 + M2*xd22 - sman*xd22 - uman*xd22 - 2*M2*xd23+ M2*xd3 + uman*xd3 + M2*xd33 - uman*xd33)))/(64*Pi**2))

! Ghost box 2
      S2NLO= cI*color*((GS**4*(tman**2*(xd1 - xd11 - 3*xd12 + 2*xd2 - 2*xd22) + tman*
     &       (4*xc0 + 2*xc2 + sman*xd0 + uman*xd0 - 2*xd00 - uman*xd1 - sman*xd11 - 2*sman*xd12 - uman*xd12 + sman*xd13 + 2*uman*xd13 - sman*xd22 - uman*xd22 + 
     &         sman*xd23 + 3*uman*xd23 - sman*xd3 - uman*xd3) + M2**2*(xd0 - xd11 - 4*xd12 + 2*xd13 + 2*xd2 - 3*xd22 + 4*xd23 - xd33) + 
     &      uman*(2*xc2 - 2*xd00 + sman*xd13 + sman*xd23 + uman*xd23 + sman*xd3 + uman*xd3 - sman*xd33 - uman*xd33) + 
     &      M2*(-4*xc0 - 4*xc2 - sman*xd0 - tman*xd0 - uman*xd0 + 4*xd00- tman*xd1 + uman*xd1 + sman*xd11 + 2*tman*xd11 + 2*sman*xd12 + 7*tman*xd12 + uman*xd12 - 
     &         2*sman*xd13 - 2*tman*xd13 - 2*uman*xd13 - 4*tman*xd2 + sman*xd22 + 5*tman*xd22 + uman*xd22 - 2*sman*xd23 - 3*tman*xd23 - 5*uman*xd23 + tman*xd3 - 
     &         uman*xd3 + sman*xd33 + 2*uman*xd33)))/(64*Pi**2))


! ---------------------  Ghost boxes u-----------------------------------
      call LTDget(M2,0d0,0d0,M2,uman,sman,M2,0d0,0d0,0d0)
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
      color = -Nc/4d0  ! only 8A

! Ghost box 1 u
      S1NLO = S1NLO + cI*color*(-(GS**4*(uman**2*(xd1 - xd11 - 3*xd12 + 2*xd2 - 2*xd22) + 
     &       uman*(4*xc0 + 2*xc2 + sman*xd0 + tman*xd0 - 2*xd00 - tman*xd1 - sman*xd11 - 2*sman*xd12 - tman*xd12 + sman*xd13 + 2*tman*xd13 - sman*xd22 - 
     &          tman*xd22 + sman*xd23 + 3*tman*xd23 - sman*xd3 - tman*xd3) + M2**2*(xd0 - xd11 - 4*xd12 + 2*xd13 + 2*xd2 - 3*xd22 + 4*xd23 - xd33) + 
     &       tman*(2*xc2 - 2*xd00 + sman*xd13 + sman*xd23 + tman*xd23 + sman*xd3 + tman*xd3 - sman*xd33 - tman*xd33) + 
     &       M2*(-4*xc0 - 4*xc2 - sman*xd0 - tman*xd0 - uman*xd0 + 4*xd00 + tman*xd1 - uman*xd1 + sman*xd11 + 2*uman*xd11 + 2*sman*xd12 + tman*xd12 + 
     &          7*uman*xd12 - 2*sman*xd13 - 2*tman*xd13 - 2*uman*xd13 - 4*uman*xd2 + sman*xd22 + tman*xd22 + 5*uman*xd22 - 2*sman*xd23 - 5*tman*xd23 - 
     &          3*uman*xd23 - tman*xd3 + uman*xd3 + sman*xd33 + 2*tman*xd33)))/(64*Pi**2))

! Ghost box 2 u
      S2NLO= S2NLO + cI*color*((GS**4*(-2*sman*xc1 - 4*M2*xc2 + M2**2*xd0 - 2*sman*xd00 - M2**2*xd11 - 2*M2**2*xd12 + M2*sman*xd12 + 2*M2**2*xd13 - M2*sman*xd2 - M2**2*xd22 + 
     &      M2*sman*xd22 + 2*M2**2*xd23 - M2*sman*xd23 + tman**2*(xd13+ xd23) + uman**2*(xd0 - xd1 + xd13 - xd2 + xd23 - xd3) - M2**2*xd33 + 
     &      uman*(-2*xc0 + 2*xc2 - 2*M2*xd0 + M2*xd1 + M2*xd11 + 2*M2*xd12 - sman*xd12 - 2*M2*xd13 + M2*xd2 + sman*xd2 + M2*xd22 - sman*xd22 - 2*M2*xd23 + 
     &         M2*xd3 + M2*xd33) + tman*(2*xc0 + 2*xc2 - M2*xd1 + uman*xd1 + M2*xd11 - uman*xd11 + 2*M2*xd12 - 2*uman*xd12 - 2*M2*xd13 - M2*xd2 + uman*xd2 + 
     &         M2*xd22 - uman*xd22 - 2*M2*xd23 + sman*xd23 - M2*xd3 + uman*xd3 + M2*xd33 - uman*xd33)))/(64*Pi**2))

! Leading order ghost amplitudes
      S1LO=(cI*GS**2*(uman -tman)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(2*sman)
      S2LO=S1LO

! Interferences between tree level and NLO
      color=Nc*CV*CF   ! given by the scalar product <c^{8A}|c^{8A}>
! S1 Interference
      M2Color(3) = 2*dble(color*dconjg(S1LO)*S1NLO)

! S2 Interference
      M2Color(3) = M2Color(3) + 2*dble(color*dconjg(S2LO)*S2NLO)  

      end 



