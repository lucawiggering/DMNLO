      double complex function DMNLO_Resum(matching,corrflag,color_basis,mi1,mi2,Wi1,Wi2,sqrtS)

      implicit none

#include "../DMNLO_Model.h"      
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

      integer nf,error_level,matching, corrflag,corrProp,corrBox,corrVert
      double precision color_Coulomb,color_basis,GS_Coulomb1,GS_Coulomb2,G_beta0,G_a1,muSc_Coulomb1,muSc_Coulomb2      
      double precision WSfav, mav, mred, mi1, mi2, Wi1, Wi2,sqrtS
      double precision mustart, muprecision
      double complex x,z,z_1,z_2,g_LO,g_LOmatch,g_NLO,g_NLOmatch1,g_NLOmatch2,G_L1,G_L2,G_v1,G_v2,G_kappa1,G_kappa2,psi_0,psi_1,psi_2
      double complex, dimension(4) :: a
      double complex, dimension(3) :: b
      double precision, parameter :: Gamma_E=0.57721566490153286d0
      double complex, parameter :: ieps=(0d0,1.d-17)                        !has to be equal to the ieps in DMNLO_LoopFunctions.F to be consistent

      
* functions      
      double complex, external :: DMNLO_DiGamma,DMNLO_TriGamma,DMNLO_PolyGamma,DMNLO_Hyper_pFq,DMNLO_HyperpFq

* functions      
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
        
        if(corrflag.eq.0) then
            corrProp = 0
            corrBox = 0
            corrVert = 0
        else
            corrProp = 999
            corrBox = 999
            corrVert = 999        
        endif


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
        call DMNLO_find_muSc(muSc_Coulomb1,GS_Coulomb1,mustart,muprecision,6,error_level,color_Coulomb,mi1,mi2,Wi1,Wi2,sqrtS)
        

        
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
     &             -2*G_kappa1*psi_0*psi_1 + 4*DMNLO_Hyper_pFq(4,3,a,b,z,1d-5,2)) + G_a1*(G_L1-psi_0 + G_kappa1*psi_1)           !Ref. [A] eq.25
     
!**************matching coefficients (for matching with full NLO- or with Coulomb O(alpha_s)-term**********
!To reduce the confusion, the g_LO above contains the LO and NLO contribution, which have to be removed eather from resummation or from NLO calculation to avoid double couting. Therefore to remove them from resummation g_NLOmatch1 and g_LOmatch were introduced. where the NLO and LO has nothing to do with the splitting into g_LO or g_NLO from above! The g_LOmatch removes the Treelevel term from the resummation, this was included in g_LO in the term with G_kappa1, the g_NLOmatch1 removes the NLO contribution from the resummation which is only included into g_LO the G_L1 term - I hope it is understandable, I was really confused by the given names

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
          g_LOmatch = (0d0,0d0)    


!use, if full Coulomb at NLO required
          g_NLOmatch1 = (0d0,0d0) 
          g_NLOmatch2 = (0d0,0d0) 

!**********************************************************************************************************

           if((corrProp.eq.999).and.(corrVert.eq.999).and.(corrBox.eq.999)) then            !NLO-Coulomb-contributions only, non-Coulomb diagrams set to zero                               

             if(matching.eq.0) then
!matching term from Coulomb O(alpha_s)-term
               DMNLO_Resum = aimag(color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*g_NLOmatch1 ) 

             else
!matching term from Coulomb O(alpha_s)-term with full NLO calculation
               DMNLO_Resum = aimag(  color_Coulomb*GS_Coulomb2**2/(4d0*pi)/G_v2*g_NLOmatch2 )  
               
             endif

           else                                                                                                            !full NLO 
 
             if(matching.eq.0) then
!matched with Coulomb O(alpha_s)-term
               
               DMNLO_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) 
     &                                          - color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*g_NLOmatch1 - g_LOmatch )     !4d0*pi/(G_v*m1**2) included to get Coulomb*LO->LO for alpha_s->0
             else if(matching.eq.2) then 

                DMNLO_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) )
             else 
!matched with full NLO calulation                                                                                                        
               DMNLO_Resum = aimag( color_Coulomb*GS_Coulomb1**2/(4d0*pi)/G_v1*(g_LO + GS_Coulomb1**2/(4d0*pi)**2*g_NLO ) 
     &                                         -  color_Coulomb*GS_Coulomb2**2/(4d0*pi)/G_v2*g_NLOmatch2 - g_LOmatch )     !4d0*pi/(G_v*m1**2) included to get Coulomb*LO->LO for alpha_s->0
             endif
             
           endif
       
             
         else

           print*,''
           print*,'Determination of muSc_Coulomb did not converge, abort evaluation!'
           print*,''

           DMNLO_Resum = 0d0

         endif


      end

!------------------------------------------------------------------------
!--- finde muSc from implicit equation and calculate GS at that scale ---          
!------------------------------------------------------------------------

      subroutine DMNLO_find_muSc(find_mu,GS_Coulomb,muSc_Coulomb,set_precision,nf,error_level,color_Coulomb,mi1,mi2,Wi1,Wi2,sqrtS)

      implicit none

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"


      integer n_count,nf,error_level
      double precision muinter,muSc_Coulomb,vel,set_precision,ASMS_nf,GS_Coulomb,color_Coulomb,color_abs
      double precision find_mu
      double precision Wsfav,mred,mav, mi1,mi2,Wi1,Wi2,sqrtS
!       double precision sqrtS

* functions      
      double precision AS,ASMS

      color_abs = dabs(color_Coulomb)   
      n_count = 0
      error_level = 0       
      muinter = muSc_Coulomb
      
      !defining average mass:
           mav = (mi1+mi2)/2d0
!defining reduced mass:
           mred = mi1*mi2/(mi1+mi2)
!defining average width
           Wsfav = (Wi1+Wi2)/2d0
           
                
!for different particles:      
      vel = dsqrt((sqrtS-2*mav)/(2d0*mred))

      if(nf.eq.6) then

      ASMS_nf = 4d0*pi*ASMS(muinter**2)       !Calculates alpha_s in the MSbar-scheme with 6 active quark flavors (including the top quark)  

      else

      ASMS_nf = 4d0*pi*AS( muinter**2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)                 !Calculates alpha_s in the MSbar-scheme with 5 active quark flavors  

      endif

!for different particles:
        do while(dabs(muinter-color_abs*ASMS_nf*2d0*mred).gt.set_precision) 

        n_count = n_count + 1

        if( muinter.gt.color_abs*ASMS_nf*2d0*mred) then

            muinter = muinter - 10d0*(muinter/(color_abs*ASMS_nf*2d0*mred)-1d0)

        else 

            muinter = muinter + 10d0*(color_abs*ASMS_nf*2d0*mred/muinter-1d0)

        endif

        if(nf.eq.6) then

        ASMS_nf = 4d0*pi*ASMS(muinter**2)                !Calculates alpha_s in the MSbar-scheme with 6 active quark flavors (including the top quark)  

        else

        ASMS_nf = 4d0*pi*AS( muinter**2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)               !Calculates alpha_s in the MSbar-scheme with 5 active quark flavors  

        endif

        
        
        if(n_count.ge.5000) then

          error_level=1
          exit

        endif

      end do


!!for different particles:
      if(muinter.lt.4d0*mred*vel) then                           !take ('more') appropriate scale to minimize higher order log(mu)-terms (Ref. [E] eq. 4.25)

        muinter = 4d0*mred*vel                                           

        if(nf.eq.6) then

          ASMS_nf = 4d0*pi*ASMS(muinter**2)            !Calculates alpha_s in the MSbar-scheme with 6 active quark flavors (including the top quark)  

        else

          ASMS_nf = 4d0*pi*AS(muinter**2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)  !Calculates alpha_s in the MSbar-scheme with 5 active quark flavors  

        endif
 
      endif
!                 

      find_mu=muinter

      GS_Coulomb = dsqrt(4*pi*ASMS_nf)                                             !Calculation of the strong coupling in the MSbar-scheme with nf=5,6 at the Coulomb-scale=muSc_Coulomb (MSbar because Greens function is calculated in MSbar, see [C] page 6) 

      end 

! ------------------------
! --- DiGamma-function ---
! ------------------------
      double complex function DMNLO_DiGamma(in)                              !from wiki(english):DiGamma-function out of 'Handbook of Mathematical Functions' p. 258 following

      implicit none


      double complex in,z,res

      z = in

      res=(0d0,0d0)

      do while(dreal(z).lt.6)

        res = res - 1d0/z

        z = z + 1d0

      end do
 
      res = res + cdlog(z) - 1d0/(2d0*z) - 1d0/(12d0*z**2) + 1d0/(120d0*z**4) - 1d0/(252d0*z**6) + 1d0/(240d0*z**8) - 5d0/(660*z**10) + 691d0/(32760d0*z**12) - 1d0/(12d0*z**14)

      DMNLO_DiGamma = res


      end function

! -------------------------
! --- TriGamma-function ---
! -------------------------
      double complex function DMNLO_TriGamma(in)                             ! out of Handbook of Mathematical Functions p. 258 following

      implicit none


      double complex in,z,res

      z = in

      res=(0d0,0d0)

      do while(dreal(z).lt.10)

        res = res + 1d0/z**2

        z = z + 1d0

      end do

      res = res + 1d0/z + 1d0/(2d0*z**2) + 1d0/(6d0*z**3) - 1d0/(30d0*z**5) + 1d0/(42d0*z**7) - 1d0/(30d0*z**9)

      DMNLO_TriGamma = res

      end function
! 
! -----------------------------
! --- PolyGamma(2)-function ---
! -----------------------------
      double complex function DMNLO_PolyGamma(in)                            !out of 'Handbook of Mathematical Functions' p. 258 following

      implicit none

      double complex in,z,res

      z = in

      res=(0d0,0d0)

      do while(dreal(z).lt.10)

        res = res - 2d0/z**3

        z = z + 1d0

      end do
 
      res = res - 1d0/z**2 - 1d0/z**3 - 1d0/(2d0*z**4) + 1d0/(6d0*z**6) - 1d0/(6d0*z**8) + 3d0/(10d0*z**10) - 5d0/(6d0*z**12) 

      DMNLO_PolyGamma = res

      end function

! ----------------------------------
! --- Logarithmic Gamma function ---
! ----------------------------------
      double complex function DMNLO_lnLanczos_Gamma(in)                      !Lancosz-Approximation of the Gamma-funtion, taken from http://en.wikipedia.org/wiki/Lanczos_approximation
                                                                       !take logarithms to be numerically more stable; while calculating the hyper-geom. function take the exponential again
      implicit none

      integer i,reflection
      double complex in,z,x,res
      integer, parameter :: g=7
      double precision, dimension(g+2) :: p
      data p/0.99999999999980993D0,676.5203681218851D0,-1259.1392167224028D0,771.32342877765313D0,-176.61502916214059D0,12.507343278686905D0,-0.13857109526572012D0,9.9843695780195716D-6,1.5056327351493116D-7/

      double precision, parameter :: pi=3.141592653589793D0

      reflection=0
      z=in

      if (dreal(z).lt.0.5d0) then 
        reflection=1
        z=1-z
      endif

      z=z-1d0

      x = p(1)       

      do i=1,g+1
        x = x+p(i+1)/(z+i)
      enddo

      res = 0.5d0*dlog(2d0*pi)+(z+0.5d0)*cdlog(z+g+0.5d0)-(z+g+0.5d0)+cdlog(x)

!Reflection formula

      if(reflection.eq.1) then

        DMNLO_lnLanczos_Gamma =  dlog(pi)-cdlog(cdsin(pi*(-z)))-res

      else

        DMNLO_lnLanczos_Gamma = res

      endif

      end function
! 
! -------------------------------------
! --- Logarithmic Pochhammer symbol ---
! -------------------------------------
      double complex function DMNLO_lnPochhammer_symbol(in,n)                !Pochhammer_symbols, taken from http://de.wikipedia.org/wiki/Pochhammer-Symbol

      implicit none

      integer n
      double complex in,res
      double complex, external :: DMNLO_lnLanczos_Gamma

        res = (0d0,0d0)

        if(n.ge.0) then

          res = DMNLO_lnLanczos_Gamma(in+n)-DMNLO_lnLanczos_Gamma(in)

        endif

        DMNLO_lnPochhammer_symbol = res

      end function
! 
! ----------------------------------------------
! --- Is the complex var a negative integer? ---
! ----------------------------------------------
      integer function DMNLO_is_negint(x)

      implicit none

      integer :: res
      double complex :: x

      double precision, parameter :: AZERO=1D-16

      if((dabs(aimag(x)).lt.azero).and.(dreal(x).lt.azero)) then

        if(dabs(idnint(dreal(x))-dreal(x)).lt.azero) then

          DMNLO_is_negint=1

          return
         
        endif

      endif

      DMNLO_is_negint=0

      end function

! ---------------------------------------------------------------------------
! --- Calculate the real part of the sum of numbers in Pochhammer-symbols ---
! ---------------------------------------------------------------------------
      double precision function DMNLO_poch_sum(p,q,a,b)

      implicit none

      integer i,p,q
      double precision res
      double complex, dimension(p) :: a
      double complex, dimension(q) :: b

      res=0d0

      do i= 1, q

        res = res + dreal(b(i))

      enddo

      do i= 1, p

        res = res - dreal(a(i))

      enddo     

      DMNLO_poch_sum = res 


      end function

! ----------------------------------------------------
! --- Hypergeometric function using series (start) ---
! ----------------------------------------------------
      double complex function DMNLO_Hyper_pFq(p,q,a,b,z,eps,improve)

      implicit none

      integer i,p,q,check
      integer, intent(in):: improve
      integer, parameter :: large_int=selected_int_kind (18)   
      integer (kind=large_int):: i_abort
      double precision eps,astop,astop_old,bstop,bstop_old
      double complex z,res

      double precision, parameter :: azero=1D-16  
      double complex, dimension(p) :: a
      double complex, dimension(q) :: b

      integer, external :: DMNLO_is_negint
      double precision, external :: DMNLO_poch_sum
      double complex, external :: DMNLO_HyperpFq,DMNLO_Hyper4F3_improve


       check=0
       i_abort=1d5

!check convergence-criteria =>see Ref. [B] appendix A-> if commented out, the convergence will be directly checked in function HyperpFq
!
!       check=1
!
!       if(p.le.q) then   
! 
!         check=0
! 
!       elseif((p.eq.q+1).and.(cdabs(z).lt.1d0)) then  
! 
!         check=0
! 
!       else if((p.eq.q+1).and.(dabs(cdabs(z)-1d0).lt.azero).and.(poch_sumstst2QQ(p,q,a,b).gt.0d0)) then   
! 
!         check=0
! 
!       else if((p.eq.q+1).and.(dabs(cdabs(z)-1d0).lt.azero).and.(dabs(dreal(z)-1d0).gt.0d0).and.(poch_sumstst2QQ(p,q,a,b).le.0d0).and.(poch_sumstst2QQ(p,q,a,b).ge.-1d0)) then  
! 
!         check=0
! 
!       endif


      if(check.eq.0) then

        do i=1,p

          if(dabs(aimag(a(i))).lt.azero) then

             if(dabs(dreal(a(i))).lt.azero) then

!               print*,'Warning in Hyper_pFq: a-Pochhammersymbol equals 0!'   

               check=2

               exit

             else if(DMNLO_is_negint(a(i)).eq.1) then

!              print*,'Warning in Hyper_pFq: a-Pochhammersymbol is a negative integer!'  

               astop=cdabs(a(i))

               if(i.eq.1) then

                 astop_old = astop

               else if(astop.lt.astop_old) then

                 astop_old = astop

               endif

             endif
         
          endif

        enddo

      endif


      if(check.eq.0) then

        do i=1,q

          if(DMNLO_is_negint(b(i)).eq.1) then

             print*,'Warning in Hyper_pFq: b-Pochhammersymbol equals 0,-1,...I will try to set a new recursion limit!'

             check=1

             bstop=cdabs(b(i))

             if(i.eq.1) then

               bstop_old = bstop

             else if(bstop.lt.bstop_old) then

               bstop_old = bstop

             endif
       
          endif

        enddo

        if((check.eq.1).and.(bstop_old.gt.astop_old)) then
          
          print*,'Warning in Hyper_pFq: new recursion-limit set!'   

          i_abort=astop_old-1 

          check=0

        endif

      endif


      if((check.eq.0).and.(improve.gt.0).and.(p.eq.4).and.(q.eq.3).and.(a(1).eq.1d0).and.(a(2).eq.1d0).and.(a(3).eq.1d0).and.(a(4).eq.1d0).and.(z.eq.1d0)) then  !improvement only written for function specified by if-statement

        res = DMNLO_Hyper4F3_improve(a,b,z,eps,improve)

      else if(check.eq.1) then
 
        print*,'Warning: Hyper_pfq-input does lead to divergent series, returning zero!'
        res = (0d0,0d0)

      else if(check.eq.2) then
 
        print*,'Warning: some Pochhammer-a equals 0, returning one (should be the right value)!'
        res = (1d0,0d0)

      else                      !standard case

        res = DMNLO_HyperpFq(p,q,a,b,z,eps,i_abort)

      endif

      DMNLO_Hyper_pFq=res

      end function
! 
! --------------------------------------------------------
! --- Hypergeometric function using series (main part) ---
! --------------------------------------------------------
      double complex function DMNLO_HyperpFq(p,q,a,b,z,eps,i_abort)          !following Ref. [B] appendix A

      implicit none

      integer i,k,l,p,q,res_fac
      double precision eps,ProgressR,ProgressIm
      double complex, dimension(p) :: a
      double complex, dimension(q) :: b
      double complex  z,res,res_old
      double complex, external :: DMNLO_lnPochhammer_symbol, DMNLO_lnLanczos_Gamma

      integer, parameter :: large_int=selected_int_kind (18)   
      integer (kind=large_int):: i_abort

      integer :: inf
      double precision :: infd
      equivalence(inf,infd)
      data inf/z'0'/     


      if(eps.lt.0d0) then 

      eps=1.d-6

      endif

      
      res_old = (0d0,0d0)  
      

      do  i=0, i_abort

      res = (0d0,0d0)   

        do k = 1,p

           res = res + DMNLO_lnPochhammer_symbol(a(k),i)

        enddo

        do k = 1,q

           res = res - DMNLO_lnPochhammer_symbol(b(k),i)

        enddo

        res = res + i*cdlog(z) - DMNLO_lnLanczos_Gamma(dcmplx(i+1))

        res  = res_old + cdexp(res)

!check progress
        progressR = dabs(dreal(res_old-res))-eps*dabs(dreal(res_old))
        progressIm = dabs(aimag(res_old-res))-eps*dabs(aimag(res_old))

        if((progressR.le.0d0).and.(progressIm.le.0d0)) exit  
  
        res_old = res

        if(i.eq.i_abort) then

!          print*,'Warning:series did not converge after',i_abort,'steps!', progressR, progressIm, eps

!          res = infd

          print*,'Warning: HyperpFq-series did not converge after',i_abort,'steps!', progressR, progressIm, eps,'returning zero!'

          res = (0d0,0d0)

        endif

      end do

      DMNLO_HyperpFq = res

      end function

! ----------------------------------------------------------
! --- Hypergeometric function using series (improvement) ---
! ----------------------------------------------------------
      double complex function DMNLO_Hyper4F3_improve(a,b,z,eps,improve)      !following Ref. [B] appendix A

      implicit none

      integer, intent(in):: improve
      double precision eps
      double complex z,res
      double precision, parameter :: azero=1D-16  
      double complex, dimension(4) :: a
      double complex, dimension(3) :: b,b1,b2,b3


      double complex, external :: DMNLO_Hyper_pFq 

     
      if(cdabs(b(1)-b(2)).lt.azero) then

        b1(1)=b(1)
        b1(2)=b(2)
        b1(3)=b(3)+1d0

        b2(1)=b(1)+1d0
        b2(2)=b(2)
        b2(3)=b(3)

        b3(1)=b(1)+1d0
        b3(2)=b(2)+1d0
        b3(3)=b(3)
      
        res = 1d0/(b(1)**2*b(3)*(b(3)-2d0*(2d0-b(1)))*(b(1)-b(3))**2) * ( b(1)**2*(b(3)-1d0)**4*DMNLO_Hyper_pFq(4,3,a,b1,z,eps,improve-1) 
     &                                        + b(1)*(b(1)-1d0)**3*b(3)*(3d0*b(1)+1d0-4d0*b(3))*DMNLO_Hyper_pFq(4,3,a,b2,z,eps,improve-1)
     &                                                         + (b(1)-1d0)**4*b(3)*(b(3)-b(1))*DMNLO_Hyper_pFq(4,3,a,b3,z,eps,improve-1)  )

      else

        b1(1)=b(1)
        b1(2)=b(2)
        b1(3)=b(3)+1d0

        b2(1)=b(1)+1d0
        b2(2)=b(2)
        b2(3)=b(3)

        b3(1)=b(1)
        b3(2)=b(2)+1d0
        b3(3)=b(3)

        res = 1d0/(b(1)+b(2)+b(3)-4d0)*( (b(1)-1d0)**4/(b(1)*(b(1)-b(2))*(b(1)-b(3)))*DMNLO_Hyper_pFq(4,3,a,b2,z,eps,improve-1)
     &                                 + (b(2)-1d0)**4/(b(2)*(b(2)-b(1))*(b(2)-b(3)))*DMNLO_Hyper_pFq(4,3,a,b3,z,eps,improve-1)
     &                                 + (b(3)-1d0)**4/(b(3)*(b(3)-b(1))*(b(3)-b(2)))*DMNLO_Hyper_pFq(4,3,a,b1,z,eps,improve-1) )

      endif

      DMNLO_Hyper4F3_improve = res

      end function DMNLO_Hyper4F3_improve
