      double complex function staustau2QQ_M2Enhancement()

      implicit none
  
#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"
#include "../util/DMNLO_RenScheme.h"
#include "../util/DMNLO_LoopIntegrals.h"

      double precision color_Coulomb,G_beta0,G_a1,muSc_Sommerfeld1,muBohr,vel    
      double complex x,z,z_1,z_2,g_LO,g_LOmatch,g_NLO,g_NLOmatch1,G_L1,G_v1,G_kappa1,psi_0,psi_1,psi_2
      double complex, dimension(4) :: a
      double complex, dimension(3) :: b
      double precision, parameter :: Gamma_E=0.57721566490153286d0
      double precision  finde_muSc
      double complex DMNLO_DiGamma,DMNLO_TriGamma,DMNLO_PolyGamma,DMNLO_Hyper_pFq

* following Ref. [A] = arXiv:0812.0919v2 and Ref. [B] = arXiv:1108.0542v2 and Ref. [C] = arXiv:0607239v1
*      and Ref. [D] = arXiv:1202.2260v2 and Ref. [E] = arXiv:1007.5414v2 and Ref. [F] = Phys. Rev. Lett.80.2531 and Ref. [G] = arXiv:1102.4295v2
                                                                            
********Couplings and scales***************

      vel = dsqrt((sqrtS-2*MSf(isf1,itt1,igen1))/MSf(isf1,itt1,igen1))
      muSc_Sommerfeld1 = 2d0*MSf(isf1,itt1,igen1)*vel
      muBohr = (EL**2/(4*pi))*MSf(isf1,itt1,igen1) ! Bohr scale

! set the Coulomb scale
      if(muBohr.gt.muSc_Sommerfeld1) then                                                  

      muSc_Sommerfeld1 = muBohr 

      endif
      
      color_Coulomb = 1d0                                              ! QED process
      G_beta0 = -4d0/3d0                                       !wikipedia: beta_function

! used in the full Greens-function
      G_v1 = cdsqrt((sqrtS+cI*WSf(isf1,itt1,igen1)-2*MSf(isf1,itt1,igen1))/MSf(isf1,itt1,igen1))  !Ref. [A] eq.26                    
      G_kappa1 = cI*EL**2/(8d0*pi*G_v1)         !Ref. [A] eq.26
      G_L1 = cdlog(cI*muSc_Sommerfeld1/(2*MSf(isf1,itt1,igen1)*G_v1))   

      psi_0 = Gamma_E + DMNLO_DiGamma(1d0-G_kappa1)   !!derivatives of psi in point 1-kappa
      psi_1 = DMNLO_TriGamma(1d0-G_kappa1)
      psi_2 = DMNLO_PolyGamma(1d0-G_kappa1)            
      
      a(1)=1d0    !POCHHAMMER SYMBOLS for GHF
      a(2)=1d0
      a(3)=1d0
      a(4)=1d0
      b(1)=2d0
      b(2)=2d0
      b(3)=1-G_kappa1
      z=1d0                                 

!************************Coefficients in Coulomb corrections***********************************************

      g_LO = -1d0/(2d0*G_kappa1)+G_L1-psi_0       
      g_NLO = 0d0
* Ref. [A] eq.25
      g_NLO =  G_beta0*(G_L1**2-2*G_L1*(psi_0-G_kappa1*psi_1) + G_kappa1*psi_2+psi_0**2-3*psi_1
     &          - 2*G_kappa1*psi_0*psi_1 + 4*DMNLO_Hyper_pFq(4,3,a,b,z,1d-5,2))+(2d0/9d0)*(G_L1-G_kappa1*psi_0+psi_1)                                    

! subtract to avoid double counting at LO (obtained by expanding the Greens function)
      g_LOmatch = cI    

      staustau2QQ_M2Enhancement = aimag( EL**2/((4d0*pi))/G_v1*(g_LO + EL**2/(4d0*pi)**2*g_NLO )  !4d0*pi/(G_v*MSf(isf1,itt1,igen1)**2) included to get Coulomb*LO->LO for alpha_s->0        
     &                                          -  g_LOmatch )    
      
      end

