      double complex function stsT2xx_Resum(matching, color_basis,mi1,mi2,Wi1,Wi2,sqrtS,colorScale)
      implicit none

#include "../../util/DMNLO_Model.h"      
#include "../../util/DMNLO_RenScheme.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer nf,error_level,matching
      double precision color_Coulomb,color_basis,GS_Coulomb1,GS_Coulomb2,G_beta0,G_a1,muSc_Coulomb1,muSc_Coulomb2      
      double precision WSfav, mav, mred, mi1, mi2, Wi1, Wi2,sqrtS,colorScale
      double precision mustart, muprecision
      double complex x,z,z_1,z_2,g_LO,g_LOmatch,g_NLO,g_NLOmatch1,g_NLOmatch2,G_L1,G_L2,G_v1,G_v2,G_kappa1,G_kappa2,psi_0,psi_1,psi_2
      double complex, dimension(4) :: a
      double complex, dimension(3) :: b
      double precision, parameter :: Gamma_E=0.57721566490153286d0
      double complex, parameter :: ieps=(0d0,1.d-17)                        !has to be equal to the ieps in DMNLO_LoopFunctions.F to be consistent


! functions      
      double complex, external :: DMNLO_DiGamma,DMNLO_TriGamma,DMNLO_PolyGamma,DMNLO_Hyper_pFq,DMNLO_HyperpFq

! functions      
      double precision AS,ASMS

!Version of Moritz, things are written for stop1stop1 incoming only!(e.g are the two incoming masses/widths taken to be the same!), for more general case see, e.g., Ref. [G]   
!Update: the incoming masses for different particles is implemted, but not checked, the calculation is implemented, but outcommented. The implementation follows the phd-thesis of Christoph Borschensky

!following Ref. [A] = arXiv:0812.0919v2
!and Ref. [B] = arXiv:1108.0542v2
!and Ref. [C] = arXiv:0607239v1
!and Ref. [D] = arXiv:1202.2260v2
!and Ref. [E] = arXiv:1007.5414v2
!and Ref. [F] = Phys. Rev. Lett.80.2531
!and Ref. [G] = arXiv:1102.4295v2


      error_level = 0

      nf=5                                                                !nf=number of quarks considered as massless at energy scale O(pcm)
                                                                      !->nf=5 standard case (see Ref. [G])


      color_Coulomb = - color_basis                 ! For the Minussign please see Diss from Christoph Borschensky P. 93 Equation (5.26) crazy definition

!defining average mass:
      mav = (mi1+mi2)/2d0
!defining reduced mass:
      mred = mi1*mi2/(mi1+mi2)
!defining average width
      Wsfav = (Wi1+Wi2)/2d0
        
!********Couplings and scales***************
      mustart = 200d0 !freely chosen 200 GeV as starting value, function solves implicite equation 
      muprecision = 0.01d0
      call DMNLO_find_muSc(muSc_Coulomb1,GS_Coulomb1,mustart, muprecision,6,error_level,colorScale,mi1,mi2,Wi1,Wi2,sqrtS)
  
                                                                      !Ref.[B] Eq.25 with accuracy 0.01 GeV for muSc and calculates 
                                                                      !gs(nf=6(->add top-quark effects to alpha_s) or nf=5) at that scale.
                                                                      !error-level says, if iteration converged. If not, abort Colomb-calculation.

      muSc_Coulomb2 = muSc                                                !Uses scale and coupling as defined in full NLO-calculation,
                                                                      !necessary for correct subtrcation of Coulomb-(O(alpha_s)) corrections from full NLO calculation
!            
      GS_Coulomb2 = GS
!******************************************

      if(error_level.eq.0) then
                        
      G_beta0 = 11d0-2d0/3d0*nf                                        !from Ref. [A] p.6 - attention CA = 3 and TF = 1/2 is used (G_beta0 = 11d0/3d0*CA-4d0/3d0*TF*nf )
      G_a1 = 31d0/3d0 - 10d0/9d0*nf                                    !from Ref. [A] p.13- attention CA = 3 and TF = 1/2 is used (G_a1 = 31d0/9d0*CA - 20d0/9d0*TF*nf)  

! !modification for different particles:  !Ref. Thesis Christoph Broschensky 

!used in the full Greens-function for different particles:
      G_v1 = cdsqrt((sqrtS+cI*WSfav-2*mav)/(2d0*mred))  
      G_kappa1 = cI*color_Coulomb*GS_Coulomb1**2/(8d0*pi*G_v1)         
      G_L1 = cdlog(cI*muSc_Coulomb1/(4*mred*G_v1))  
     
!used for matching the loop-part in NLO-calculation for different particles:
      G_v2 = cdsqrt((sqrtS+ieps-2*mav)/(2d0*mred))                              ! -> chose width=ieps as used in DMNLO_Loop.F
      G_kappa2 = cI*color_Coulomb*GS_Coulomb2**2/(8d0*pi*G_v2)         
      G_L2 = cdlog(cI*muSc_Coulomb2/(4*mred*G_v2))  

      psi_0 = Gamma_E + DMNLO_DiGamma(1d0-G_kappa1) 
      psi_1 = DMNLO_TriGamma(1d0-G_kappa1)
      psi_2 = DMNLO_PolyGamma(1d0-G_kappa1)            
 
!--------------Set Arguments of Hypergeometric function----------------------------------------------       
      a(1)=1d0
      a(2)=1d0
      a(3)=1d0
      a(4)=1d0

      b(1)=2d0
      b(2)=2d0
      b(3)=1d0-G_kappa1

      z=1d0                                 

!************************Coefficients in Coulomb corrections***********************************************

      g_LO = -1d0/(2d0*G_kappa1)+G_L1-psi_0                                                                           !Ref. [A] eq.25

      g_NLO =  G_beta0*(G_L1**2-2*G_L1*(psi_0-G_kappa1*psi_1) + G_kappa1*psi_2+psi_0**2-3*psi_1
     &            -2*G_kappa1*psi_0*psi_1 + 4*DMNLO_Hyper_pFq(4,3,a,b,z,1d-5,2)) + G_a1*(G_L1-psi_0 + G_kappa1*psi_1)           !Ref. [A] eq.25

!**************matching coefficients (for matching with full NLO- or with Coulomb O(alpha_s)-term**********
!To reduce the confusion, the g_LO above contains the LO and NLO contribution, which have to be removed eather from resummation 
!or from NLO calculation to avoid double couting. Therefore to remove them from resummation g_NLOmatch1 and g_LOmatch were introduced,
! where the NLO and LO has nothing to do with the splitting into g_LO or g_NLO from above!
!The g_LOmatch removes the Treelevel term from the resummation, this was included in g_LO in the term with G_kappa1, the g_NLOmatch1 removes 
!the NLO contribution from the resummation which is only included into g_LO the G_L1 term - I hope it is understandable, I was really confused by the given names

! subtract to avoid double counting at LO 
      g_LOmatch = cI    
!for different particles:
!subtract from Greens-function: Ref. [A] eq.27, matching-coefficient to avoid double counting at NLO           
      g_NLOmatch1 = (cI*pi/2d0 + dlog(muSc_Coulomb1/(4*mred)) - cdlog(G_v1))   

!subtract from full NLO: Ref. [A] eq.27, matching-coefficient to avoid double counting at NLO
      g_NLOmatch2 = (cI*pi/2d0 + dlog(muSc_Coulomb2/(4*mred)) - cdlog(G_v2))   
!          

***************matching coefficient if full Coulomb correction is required*********************************

! !use, if full Coulomb at LO required
!    g_LOmatch = (0d0,0d0)    


!use, if full Coulomb at NLO required
!    g_NLOmatch1 = (0d0,0d0) 
!    g_NLOmatch2 = (0d0,0d0) 
!          print*,'corrProp etc',corrBox,corrProp,corrVert
!**********************************************************************************************************

      if(matching.eq.0) then
!matched with Coulomb O(alpha_s)-term
         
            stsT2xx_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) 
     &                                          - color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*g_NLOmatch1 - g_LOmatch )     !4d0*pi/(G_v*m1**2) included to get Coulomb*LO->LO for alpha_s->0

      else if(matching.eq.2) then
            stsT2xx_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) )
      else                                       
!matched with full NLO calulation                                                                                                        
            stsT2xx_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) 
     &                                         -  color_Coulomb*GS_Coulomb2**2/(4d0*pi)/G_v2*g_NLOmatch2 - g_LOmatch )     !4d0*pi/(G_v*m1**2) included to get Coulomb*LO->LO for alpha_s->0
      endif
 
       
      else

      print*,''
      print*,'Determination of muSc_Coulomb did not converge, abort evaluation!'
      print*,''

      stsT2xx_Resum = 0d0

      endif



      end

