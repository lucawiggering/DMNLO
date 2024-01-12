!**********************************************************************************************************************************

      subroutine QQ2xx_M2onshellsubtract(argx,argy,argxh,argyh,Weight)

      implicit none

#include "QQ2xx_Model.h"
#include "QQ2xx_Kinematics.h"
	
      integer argx,argy,argxh,argyh
      double precision Weight,Width_elw,Width_tot

      external QQ2xx_widthS, QQ2xx_widthV

      Weight = 0d0
      Width_elw = 0d0
      Width_tot = 0d0
 
      if((abort22.eq.2).and.((argx.eq.1).or.(argy.eq.1))) then      !there must be one Gluon in the final state

        if((argx.eq.1).and.(argyh.ne.0)) then                       !Gluon-Scalar (Scalar-decay)

          call QQ2xx_widthS(argyh,Width_elw,Width_tot)

        else if((argxh.ne.0).and.(argy.eq.1)) then                  !Scalar-Gluon (Scalar-decay)

          call QQ2xx_widthS(argxh,Width_elw,Width_tot)

        else if((argx.eq.1).and.(argy.eq.3)) then                   !Gluon-Z (Z-decay, charged final state not yet covered (e.g. W+ G))

          call QQ2xx_widthV(argy,Width_elw,Width_tot)

        else if((argx.eq.3).and.(argy.eq.1)) then                   !Z-Gluon (Z-decay, charged final state not yet covered (e.g. W+ G))

          call QQ2xx_widthV(argx,Width_elw,Width_tot)

        endif
 
        if(dabs(Width_tot).gt.zerocheck) then

           Weight = Width_elw/Width_tot

        else

           Weight = 0d0                                             !particle is stable -> no subtration necessary

        endif

      endif

!      print*,'final Weight',argx,argy,argxh,argyh,Weight
                                                                                                                              
!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_widthS(i,WidthS_elw,WidthS_tot)              !calculates width for scalar decay into elw. and into all possible final states (only decay of neutral scalars included so far, since charged final states are not yet covered)

      implicit none

#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Init2to2.h"

      integer i,argx,argy,gen
      double precision s_in,m_a,m_b,int_factor,Symm_factorS,Symm_factorV,Width_res,WidthS_elw,WidthS_tot,MHVV_width,
     &                 MHGhGh_width,MHVS_width,MHSV_width,MHSS_width,MHFFbar_width
      complex*16 Coupl_1,Coupl_2,Coupl_1L,Coupl_2L,Coupl_1R,Coupl_2R


      s_in = MHk(i)**2                    !ingoing energy = mass of decaying particle 
      WidthS_elw = 0d0
      WidthS_tot = 0d0

      do argx=3,4                         !no Gluon/Photon-Higgs coupling -> start at argx=3 (photon-coupling to Higgs only for H+/- -> but charged scalar decay not yet included)
         argy=argx                        !avoid double counting, no charged final state included yet
!        do argy=argx,4                    !avoid double counting of final states

!Vector-Vector final states     
             
          m_a = MV(argx)                  !outgoing mass a (->k1)
          m_b = MV(argy)                  !outgoing mass b (->k2)
          Coupl_1 = LLR(argx,argy,i)
          Coupl_2 = dconjg(LLR(argx,argy,i))    

          Width_res = int_factor(s_in,m_a,m_b)*Symm_factorV(argx,argy)*MHVV_width(s_in,m_a,m_b,Coupl_1,Coupl_2)          
          
          WidthS_elw =  WidthS_elw + Width_res
          WidthS_tot =  WidthS_tot + Width_res

!Ghost-Ghost final states

          m_a = MV(argx) 
          m_b = MV(argy)                        
          Coupl_1 = GGH(argx,argy,i)
          Coupl_2 = GGH(argy,argx,i) 
          
          Width_res = int_factor(s_in,m_a,m_b)*Symm_factorV(argx,argy)*MHGhGh_width(s_in,m_a,m_b,Coupl_1,Coupl_2)

	  WidthS_elw = WidthS_elw + Width_res
	  WidthS_tot = WidthS_tot + Width_res

!        enddo  

        do argy=1,6 

!Vector-Scalar final states
                     
          m_a = MV(argx) 
          m_b = MHk(argy) 
          Coupl_1 = JLR(argx,i,argy)
          Coupl_2 = dconjg(JLR(argx,i,argy))  

          if((argx.eq.4).and.(argy.eq.6))  then                                                       

            Width_res = 2*int_factor(s_in,m_a,m_b)*MHVS_width(s_in,m_a,m_b,Coupl_1,Coupl_2)    !add up Goldstone boson contributions for  W+W- final state
          
            WidthS_elw = WidthS_elw +  Width_res
            WidthS_tot = WidthS_tot +  Width_res

          else
          
            Width_res = int_factor(s_in,m_a,m_b)*MHVS_width(s_in,m_a,m_b,Coupl_1,Coupl_2)      !double counting of ZZ final state canceled against symmetry factor

            WidthS_elw = WidthS_elw + Width_res
            WidthS_tot = WidthS_tot + Width_res

          endif 

        enddo
      enddo

!Scalar-Scalar final states

      do argx=1,6
         do argy=argx,6
        
          m_a = MHk(argx)         
          m_b = MHk(argy) 
          Coupl_1 = ILR(i,argx,argy)
          Coupl_2 = dconjg(ILR(i,argx,argy))
          
          Width_res = int_factor(s_in,m_a,m_b)*Symm_factorS(argx,argy)*MHSS_width(s_in,m_a,m_b,Coupl_1,Coupl_2)

          WidthS_elw = WidthS_elw + Width_res
          WidthS_tot = WidthS_tot + Width_res 

        enddo
      enddo
 
!Fermion-Antifermion final states

      do argx=1,4    
        do argy=1,4              
          do gen=1,3                      !no generation mixing

              m_a = Mf(argx,gen) 
              m_b = Mf(argy,gen)
              Coupl_1L = DL(i,argy,argx,gen,gen)
              Coupl_1R = DR(i,argy,argx,gen,gen)
              Coupl_2L = dconjg(DL(i,argy,argx,gen,gen))
              Coupl_2R = dconjg(DR(i,argy,argx,gen,gen))

              if((argx.le.2).and.(argy.le.2)) then                  !no mixing between elw. and strong final states -> mixed final states not included
              
                Width_res = int_factor(s_in,m_a,m_b)*MHFFbar_width(s_in,m_a,m_b,Coupl_1L,Coupl_1R,Coupl_2L,Coupl_2R)

                WidthS_elw = WidthS_elw + Width_res
                WidthS_tot = WidthS_tot + Width_res

              else if ((argx.gt.2).and.(argy.gt.2)) then
              
                Width_res = int_factor(s_in,m_a,m_b)*NC*MHFFbar_width(s_in,m_a,m_b,Coupl_1L,Coupl_1R,Coupl_2L,Coupl_2R) 

                WidthS_tot = WidthS_tot + Width_res

              endif

          enddo
        enddo
      enddo

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_widthV(i,WidthV_elw,WidthV_tot)              !calculates width for vector decay into elw. and into all possible final states (only decay of neutral vectors included so far)

      implicit none

#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Init2to2.h"

      integer i,argx,argy,gen
      double precision s_in,m_a,m_b,int_factor,Width_res,WidthV_elw,WidthV_tot,MVFFbar_width
      complex*16 Coupl_1,Coupl_2,Coupl_1L,Coupl_2L,Coupl_1R,Coupl_2R


      s_in = MV(i)**2                     !ingoing energy = mass of decaying particle -> only fermionic final states accesable for on-shell-decay
      WidthV_elw = 0d0
      WidthV_tot = 0d0

!Fermion-Antifermion final states

      do argx=1,4    
        do argy=1,4              
          do gen=1,3                      !no generation mixing

              m_a = Mf(argx,gen) 
              m_b = Mf(argy,gen)
              Coupl_1L = BL(i,argy,argx,gen,gen)
              Coupl_1R = BR(i,argy,argx,gen,gen)
              Coupl_2L = dconjg(BL(i,argy,argx,gen,gen))
              Coupl_2R = dconjg(BR(i,argy,argx,gen,gen))


              if((argx.le.2).and.(argy.le.2)) then                  !no mixing between elw. and strong final states -> mixed final states not included

                Width_res = int_factor(s_in,m_a,m_b)*MVFFbar_width(s_in,m_a,m_b,Coupl_1L,Coupl_1R,Coupl_2L,Coupl_2R)
              
                WidthV_elw = WidthV_elw + Width_res
                WidthV_tot = WidthV_tot + Width_res

              else if ((argx.gt.2).and.(argy.gt.2)) then
              
                Width_res = int_factor(s_in,m_a,m_b)*NC*MVFFbar_width(s_in,m_a,m_b,Coupl_1L,Coupl_1R,Coupl_2L,Coupl_2R)

                WidthV_tot = WidthV_tot + Width_res 

              endif

          enddo
        enddo
      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      double precision function Symm_factorS(argx,argy)

      integer argx,argy
      double precision factor


      if((argx.ne.5).and.(argx.ne.6).and.(argx.eq.argy)) then

        factor=0.5d0

      else

        factor=1d0
 
      endif
 
      Symm_factorS = factor

!==================================================================================================================================

      end function
    

!**********************************************************************************************************************************
!**********************************************************************************************************************************


      double precision function Symm_factorV(argx,argy)

      integer argx,argy
      double precision factor


      if((argx.ne.4).and.(argx.eq.argy)) then

        factor=0.5d0

      else

        factor=1d0

      endif

      Symm_factorV = factor

!==================================================================================================================================

      end function
      

!**********************************************************************************************************************************
!**********************************************************************************************************************************


      double precision function int_factor(s_in,m_a,m_b)

#include "QQ2xx_Model.h"
#include "QQ2xx_Kinematics.h"

      double precision s_in,m_a,m_b,kappa,factor

      if(dsqrt(s_in).ge.m_a+m_b) then

        factor = kappa(s_in,m_a**2,m_b**2)/s_in

      else

        factor = 0d0                      !set intfactor to zero if final state is kinematially not accessable

      endif

      int_factor = factor

!==================================================================================================================================

      end function

      