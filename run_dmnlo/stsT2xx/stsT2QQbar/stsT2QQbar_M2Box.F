! Results for the boxes produced with FormCalc 9.9. These are in perfect agreement with the results below.

!      double precision function FormCalcBox1()
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call LTDget(sman,0d0,tman,m1**2,0d0,m1**2,0d0,0d0,0d0,m1**2) 
!      call myCget(sman,0d0,0d0,0d0,0d0,0d0)
!
!      FormCalcBox1 = 2*dble((-7*gs**6*(m1**4 - tman*uman)*(4*xc0 - 4*tman*xd0 - 2*xd00 + tman*xd23 - 4*tman*xd3 + m1**2*(4*xd0 - xd23 - xd33) - tman*xd33))/(12.*Pi**2*sman))     
!
!      end      
!
!      double precision function FormCalcBox1u()
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call LTDget(sman,0d0,uman,m1**2,0d0,m1**2,0d0,0d0,0d0,m1**2)
!      call myCget(sman,0d0,0d0,0d0,0d0,0d0)
!
!      FormCalcBox1u = 2*dble((gs**6*(-m1**4 + tman*uman)*((-2*(-m1**2 + uman)*xd0)/3. - xd00/3. + (2*(xc0 - uman*xd3))/3. + ((-m1**2 + uman)*xd23 + (-m1**2 - uman)*xd33)/6.))/(Pi**2*sman))     
!
!      end    
!
!      double precision function FormCalcBox2(fType,fGen)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      integer fType,fGen,k
!
!      FormCalcBox2 = 0d0
!
!      do k=1,2
!
!      call LTDget(sman,m1**2,tman,0d0,m1**2,0d0,Mgl2,Mgl2,Mf2(itt1,igen1),MSf2(k,fType,fGen))
!      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))
!
!      FormCalcBox2 = FormCalcBox2 + 2*dble( (7*gs**6*(-m1**4 + tman*uman)*(xc0*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     &         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)) + 
!     &      (-m1**2 + tman)*xd3*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     &         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)) - 
!     &      Mgl**2*xd0*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     &         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen) + 
!     &         ((m1**2 - Mf2(itt1,igen1))*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     &             USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2) + 
!     &      xd2*(Mgl*Mf(itt1,igen1)*USf(isf1,1,itt1,igen1)*USfC(isf1,2,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
!     &         Mgl*(Mf(itt1,igen1)*USf(isf1,2,itt1,igen1)*USfC(isf1,1,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
!     &            Mgl*(-(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen)) - 
!     &               USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen) - 
!     &               ((m1**2 + Mf2(itt1,igen1))*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     &                    USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2)))))/(6.*Pi**2*sman))     
!
!      enddo
!
!      end    
!
!      double precision function FormCalcBox2u(fType,fGen)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      integer fType,fGen,k
!
!      FormCalcBox2u = 0d0
!
!      do k=1,2
!
!      call LTDget(sman,m1**2,uman,0d0,m1**2,0d0,Mgl2,Mgl2,Mf2(itt1,igen1),MSf2(k,fType,fGen))
!      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))
!
!      FormCalcBox2u = FormCalcBox2u + 2*dble((gs**6*(-m1**4 + tman*uman)*(xc0*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     -         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)) + 
!     -      (-m1**2 + uman)*xd3*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     -         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)) - 
!     -      Mgl**2*xd0*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     -         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen) + 
!     -         ((m1**2 - Mf2(itt1,igen1))*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     -             USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2) + 
!     -      xd2*(Mgl*Mf(itt1,igen1)*USf(isf1,1,itt1,igen1)*USfC(isf1,2,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
!     -         Mgl*(Mf(itt1,igen1)*USf(isf1,2,itt1,igen1)*USfC(isf1,1,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
!     -            Mgl*(-(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen)) - 
!     -               USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen) - 
!     -               ((m1**2 + Mf2(itt1,igen1))*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
!     -                    USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2)))))/(3.*Pi**2*sman))     
!
!      enddo
!
!      end 
      
      double precision function stsT2QQbar_M2Box()
	implicit none
! This function gives the contribution of the octet part to the squared matrix element. 
! The singlet part of the squared matrix element is always zero since the tree level (s-channel) amplitude lives only in the octet subspace.
! The variable "colori" is the coefficient of \delta_{si} \delta_{tj}. (The part \delta_{st} \delta_{ij}  of the octet colour basis element
! vanishes under contraction with the tree amplitude.)

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_LoopIntegrals.h"

	integer k,fType,fGen
	double complex  g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3,g3L,g3R,DenSV
	double precision Krondelta,M0l,M1l,M2l,M3l
	double complex vbA, vbBL, vbBR, vbCLp1, vbCRp1, vbCLk1, vbCRk1, vbCLk2, vbCRk2, vbDp1, vbDk1, vbDk2
	double precision colorf, colori

!      double precision FormCalcBox1u,FormCalcBox1,FormCalcBox2,FormCalcBox2u


! Gluon propagator
	DenSV = sqrtS**2

! initilize with zero
      stsT2QQbar_M2Box = 0d0

! Loop over massless quarks
      do fType = 3,4
      do fGen = 1,2

! Set all form factors to zero.	
      vbA = 0d0
      vbBL = 0d0
      vbBR = 0d0
      vbCLp1 = 0d0
      vbCRp1 = 0d0
      vbCLk1 = 0d0
      vbCRk1 = 0d0
      vbCLk2 = 0d0
      vbCRk2 = 0d0
      vbDp1 = 0d0
      vbDk1 = 0d0
      vbDk2 = 0d0

	
! ***************************************************************************
! ************ Box 1: Squark-Gluon-Quark-Gluon  
! ***************************************************************************
      g0 = BL(fType,fType,fGen,fGen)
      g1 = SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
      g2 = SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
      g3 = BL(fType,fType,fGen,fGen)
   
      call LTDget(m1**2,0d0,0d0,m1**2,tman,sman, m1**2,0d0,0d0,0d0)
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
		  		  
	colori =(Nc**2-2d0)/(4d0*Nc)

	vbDp1 = vbDp1 + colori*g0*g1*g2*g3*(4*xc0 - tman*xd0 - 2*xd00 + 2*tman*xd1 - tman*xd11 - 3*tman*xd12 - sman*xd13 - 3*tman*xd13 - uman*xd13 + tman*xd2 - 2*tman*xd22 - 
     -    sman*xd23 - 4*tman*xd23 - uman*xd23 - sman*xd3 + tman*xd3 - uman*xd3 - sman*xd33 - 2*tman*xd33 - uman*xd33 + 
     -    m1**2*(3*xd0 - 2*xd1 - xd11 - xd12 - xd2 + xd23 + xd33))


! ***************************************************************************
! ************ Box 1 u: Squark-Gluon-Quark-Gluon ("u-channel")
! ***************************************************************************

      g0 = BL(fType,fType,fGen,fGen)
      g1 = SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
      g2 = SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
      g3 = BL(fType,fType,fGen,fGen)

      call LTDget(m1**2,0d0,0d0,m1**2,uman,sman, m1**2,0d0,0d0,0d0)
      call myCget(0d0,0d0,sman,0d0,0d0,0d0)
                        
      colori = - 1d0/(2*Nc)
	
	vbDp1 = vbDp1 + colori*g0*g1*g2*g3*(-4*xc0 + uman*xd0 + 2*xd00 - 2*uman*xd1 + uman*xd11 + 3*uman*xd12 + sman*xd13 + tman*xd13 + 3*uman*xd13 - uman*xd2 + 
     -    2*uman*xd22 + sman*xd23 + tman*xd23 + 4*uman*xd23 + sman*xd3 + tman*xd3 - uman*xd3 + 
     -    m1**2*(-3*xd0 + 2*xd1 + xd11 + xd12 + xd2 - xd23 - xd33) + sman*xd33 + tman*xd33 + 2*uman*xd33)

	
! ***************************************************************************
! ************ Box 2: Quark-Gluino-Squark(k)-Gluino 
! ***************************************************************************
 
      do k = 1,2 ! Sfermion index
		
      M0l = Mf(itt1,igen1)
      M1l = MGl
      M2l = MSf(k,fType,fGen)
      M3l = MGl

      g0L = GluL(k,fType,fType,fGen,fGen)
      g0R = GluR(k,fType,fType,fGen,fGen)
      g1L = GluL(isf2,itt2,itt1,igen2,igen1)
      g1R = GluR(isf2,itt2,itt1,igen2,igen1)
      g2L = GluR(isf1,itt2,itt1,igen2,igen1)
      g2R = GluL(isf1,itt2,itt1,igen2,igen1)
      g3L = GluR(k,fType,fType,fGen,fGen)
      g3R = GluL(k,fType,fType,fGen,fGen)

      call LTDget(m1**2,0d0,0d0,m1**2,tman,sman, M0l**2,M1l**2,M2l**2,M3l**2)
      call myCget(0d0,0d0,sman,M1l**2,M2l**2,M3l**2)
		  		  
	colori = (Nc**2-2d0)/(4*Nc)
	
	vbBL = vbBL + colori*(g0L*g3L*(g1L*g2L*M0l*M1l*M3l*xd0 + g1L*g2R*M3l*(xc0 + M0l**2*xd0 + m1**2*xd1 + tman*xd2 + tman*xd3) + 
     -    g1R*g2L*M1l*(xc0 + M0l**2*xd0 + tman*xd1 + tman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3) + 
     -    g1R*g2R*M0l*(xc0 + M0l**2*xd0 + tman*xd0 + m1**2*xd1 + tman*xd1 + 2*tman*xd2 - m1**2*xd3 + sman*xd3 + 2*tman*xd3 + uman*xd3)))
	
	vbBR = vbBR + colori*(g0R*g3R*(g1R*g2R*M0l*M1l*M3l*xd0 + g1R*g2L*M3l*(xc0 + M0l**2*xd0 + m1**2*xd1 + tman*xd2 + tman*xd3) + 
     -    g1L*g2R*M1l*(xc0 + M0l**2*xd0 + tman*xd1 + tman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3) + 
     -    g1L*g2L*M0l*(xc0 + M0l**2*xd0 + tman*xd0 + m1**2*xd1 + tman*xd1 + 2*tman*xd2 - m1**2*xd3 + sman*xd3 + 2*tman*xd3 + uman*xd3)))
     
	vbCLp1 = vbCLp1 + colori*(g0L*g3R*(g1L*g2R*M1l*M3l*(xd1 + xd2 + xd3) + g1R*g2R*M0l*M1l*(xd0 + xd1 + xd2 + xd3) + g1L*g2L*M0l*M3l*(xd0 + xd1 + xd2 + xd3) + 
     -    g1R*g2L*(xc0 + m1**2*xd1 + tman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 + M0l**2*(2*xd0 + xd1 + xd2 + xd3))))
	
	vbCRp1 = vbCRp1 + colori*(g0R*g3L*(g1R*g2L*M1l*M3l*(xd1 + xd2 + xd3) + g1L*g2L*M0l*M1l*(xd0 + xd1 + xd2 + xd3) + g1R*g2R*M0l*M3l*(xd0 + xd1 + xd2 + xd3) + 
     -    g1L*g2R*(xc0 + m1**2*xd1 + tman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 + M0l**2*(2*xd0 + xd1 + xd2 + xd3))))
		
      enddo


! ***************************************************************************
! ************ Box 2 u: Quark-Gluino-Squark(k)-Gluino ("u-channel")
! ***************************************************************************

      do k = 1,2 ! Sfermion index

      M0l = Mf(itt1,igen1)
      M1l = MGl
      M2l = MSf(k,fType,fGen)
      M3l = MGl

      g0L = GluR(k,fType,fType,fGen,fGen)
      g0R = GluL(k,fType,fType,fGen,fGen)
      g1L = GluL(isf2,itt2,itt1,igen2,igen1)
      g1R = GluR(isf2,itt2,itt1,igen2,igen1)
      g2L = GluR(isf1,itt2,itt1,igen2,igen1)
      g2R = GluL(isf1,itt2,itt1,igen2,igen1)
      g3L = GluL(k,fType,fType,fGen,fGen)
      g3R = GluR(k,fType,fType,fGen,fGen)

      call LTDget(m1**2,0d0,0d0,m1**2,uman,sman, M0l**2,M1l**2,M2l**2,M3l**2)
	call myCget(0d0,0d0,sman,M1l**2,M2l**2,M3l**2)

      colori = - 1d0/(2*Nc)

      vbBL = vbBL + colori*(g0L*g3L*(g1L*g2L*M0l*M1l*M3l*xd0 + g1L*g2R*M3l*(xc0 + M0l**2*xd0 + m1**2*xd1 + uman*xd2 + uman*xd3) + 
     -    g1R*g2L*M1l*(xc0 + M0l**2*xd0 + uman*xd1 + uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3) + 
     -    g1R*g2R*M0l*(xc0 + M0l**2*xd0 + uman*xd0 + m1**2*xd1 + uman*xd1 + 2*uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + 2*uman*xd3)))

      vbBR = vbBR + colori*(g0R*g3R*(g1R*g2R*M0l*M1l*M3l*xd0 + g1R*g2L*M3l*(xc0 + M0l**2*xd0 + m1**2*xd1 + uman*xd2 + uman*xd3) + 
     -    g1L*g2R*M1l*(xc0 + M0l**2*xd0 + uman*xd1 + uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3) + 
     -    g1L*g2L*M0l*(xc0 + M0l**2*xd0 + uman*xd0 + m1**2*xd1 + uman*xd1 + 2*uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + 2*uman*xd3)))

      vbCLp1 = vbCLp1 + colori*(-(g0R*g3L*(g1R*g2L*M1l*M3l*(xd1 + xd2 + xd3) + g1L*g2L*M0l*M1l*(xd0 + xd1 + xd2 + xd3) + g1R*g2R*M0l*M3l*(xd0 + xd1 + xd2 + xd3) + 
     -     g1L*g2R*(xc0 + m1**2*xd1 + uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 + M0l**2*(2*xd0 + xd1 + xd2 + xd3)))))

      vbCRp1 = vbCRp1 + colori*(-(g0L*g3R*(g1L*g2R*M1l*M3l*(xd1 + xd2 + xd3) + g1R*g2R*M0l*M1l*(xd0 + xd1 + xd2 + xd3) + g1L*g2L*M0l*M3l*(xd0 + xd1 + xd2 + xd3) + 
     -     g1R*g2L*(xc0 + m1**2*xd1 + uman*xd2 - m1**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 + M0l**2*(2*xd0 + xd1 + xd2 + xd3)))))

      enddo
  

! -------------------------------------------------------------------------------------------------
! -----------------here take /(4PI^2) with ---------------------------------------------------------
! -------------------------------------------------------------------------------------------------
	  
      vbA = vbA/(4d0*Pi)**2
      vbBL = vbBL/(4d0*Pi)**2
      vbBR = vbBR/(4d0*Pi)**2
      vbCLp1 = vbCLp1/(4d0*Pi)**2
      vbCRp1 = vbCRp1/(4d0*Pi)**2
      vbCLk1 = vbCLk1/(4d0*Pi)**2
      vbCRk1 = vbCRk1/(4d0*Pi)**2
      vbCLk2 = vbCLk2/(4d0*Pi)**2
      vbCRk2 = vbCRk2/(4d0*Pi)**2
      vbDp1 = vbDp1/(4d0*Pi)**2
      vbDk1 = vbDk1/(4d0*Pi)**2
      vbDk2 = vbDk2/(4d0*Pi)**2

! c--------------------------------------------------------------------------------------------------------------------------
! c----------------------------------------------amplitude squared----------------------------------------------------------
! c--------------------------------------------------------------------------------------------------------------------------


!----------------------------------------------Box-S-channel Gluon----------------------------------------------------------
      colorf = 1d0/2d0*(Nc**2-1)

      stsT2QQbar_M2Box =  stsT2QQbar_M2Box - 2d0*colorf*dreal((sqrtS**3*((muu4**3*(vbA + vbBR) + muu4*(-1 + 2*muu1**2 - 2*tred)*(vbA + vbBR) + muu3**3*(vbA + vbBL + muu4*sqrtS*(vbCRk1 - vbCRk2 + vbDk1 - vbDk2)) - 
     -           muu4**4*sqrtS*(vbCLk2 + vbDk2) + 2*sqrtS*(muu1**2*(muu2**2 - tred) + tred*(1 - muu2**2 + tred))*(vbCLp1 + vbDp1) + muu3**4*sqrtS*(vbCLk1 + vbCLp1 + vbDk1 + vbDp1) - 
     -           muu4**2*sqrtS*((-1 + 2*muu1**2 - 2*tred)*vbCLk2 + tred*vbCLp1 - vbDk2 - 2*tred*vbDk2 + tred*vbDp1 + muu1**2*(vbCLp1 + 2*vbDk2 + vbDp1)) + 
     -           muu3*((-1 + 2*muu2**2 + muu4**2 - 2*tred)*vbA + (-1 + 2*muu2**2 + muu4**2 - 2*tred)*vbBL + 
     -              muu4*sqrtS*((-1 + 2*muu1**2 + muu4**2 - 2*tred)*vbCRk1 + vbCRk2 - 2*muu2**2*vbCRk2 - muu4**2*vbCRk2 + 2*tred*vbCRk2 - vbCRp1 + 3*muu1**2*vbCRp1 + muu2**2*vbCRp1 - 
     -                 vbDk1 + 2*muu1**2*vbDk1 + muu4**2*vbDk1 - 2*tred*vbDk1 + vbDk2 - 2*muu2**2*vbDk2 - muu4**2*vbDk2 + 2*tred*vbDk2 - vbDp1 + 3*muu1**2*vbDp1 + muu2**2*vbDp1)) + 
     -           muu3**2*(muu4*(vbA + vbBR) + muu4**2*sqrtS*(vbCLk1 - vbCLk2 + vbCLp1 + vbDk1 - vbDk2 + vbDp1) + 
     -              sqrtS*((-1 + 2*muu2**2 - 2*tred)*vbCLk1 + (-1 + muu2**2 - 3*tred)*vbCLp1 - vbDk1 + 2*muu2**2*vbDk1 - 2*tred*vbDk1 - vbDp1 + muu2**2*vbDp1 - 3*tred*vbDp1)))*
     -         BL(fType,fType,fGen,fGen) + (muu4**3*(vbA + vbBL) + muu4*(-1 + 2*muu1**2 - 2*tred)*(vbA + vbBL) + muu3**3*(vbA + vbBR + muu4*sqrtS*(vbCLk1 - vbCLk2 + vbDk1 - vbDk2)) - 
     -           muu4**4*sqrtS*(vbCRk2 + vbDk2) + 2*sqrtS*(muu1**2*(muu2**2 - tred) + tred*(1 - muu2**2 + tred))*(vbCRp1 + vbDp1) + muu3**4*sqrtS*(vbCRk1 + vbCRp1 + vbDk1 + vbDp1) - 
     -           muu4**2*sqrtS*((-1 + 2*muu1**2 - 2*tred)*vbCRk2 + tred*vbCRp1 - vbDk2 - 2*tred*vbDk2 + tred*vbDp1 + muu1**2*(vbCRp1 + 2*vbDk2 + vbDp1)) + 
     -           muu3*((-1 + 2*muu2**2 + muu4**2 - 2*tred)*vbA + (-1 + 2*muu2**2 + muu4**2 - 2*tred)*vbBR + 
     -              muu4*sqrtS*((-1 + 2*muu1**2 + muu4**2 - 2*tred)*vbCLk1 + vbCLk2 - 2*muu2**2*vbCLk2 - muu4**2*vbCLk2 + 2*tred*vbCLk2 - vbCLp1 + 3*muu1**2*vbCLp1 + muu2**2*vbCLp1 - 
     -                 vbDk1 + 2*muu1**2*vbDk1 + muu4**2*vbDk1 - 2*tred*vbDk1 + vbDk2 - 2*muu2**2*vbDk2 - muu4**2*vbDk2 + 2*tred*vbDk2 - vbDp1 + 3*muu1**2*vbDp1 + muu2**2*vbDp1)) + 
     -           muu3**2*(muu4*(vbA + vbBL) + muu4**2*sqrtS*(vbCRk1 - vbCRk2 + vbCRp1 + vbDk1 - vbDk2 + vbDp1) + 
     -              sqrtS*((-1 + 2*muu2**2 - 2*tred)*vbCRk1 + (-1 + muu2**2 - 3*tred)*vbCRp1 - vbDk1 + 2*muu2**2*vbDk1 - 2*tred*vbDk1 - vbDp1 + muu2**2*vbDp1 - 3*tred*vbDp1)))*
     -         BR(fType,fType,fGen,fGen))*dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/dconjg(DenSV))


      enddo
      enddo

      end


