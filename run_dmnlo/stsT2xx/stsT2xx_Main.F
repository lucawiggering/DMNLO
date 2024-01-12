* DM@NLO Project - B.Herrmann, K.Kovarik, M.Klasen, Q.Le Boulc`h, J.Harz, P.Steppeler, M. Meinecke, S. Schmiemann, J.Branahl, L. Wiggering

      subroutine stst2xx(renscale,pcm,iflag,myresult,myerror,gauge,gauge2,slhafilename)
      implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_GenCouplings.h" 
#include "common/stsT2xx_GenCounterterms.h"
#include "../util/DMNLO_RenScheme.h"

	integer ndim, ncomp, iflag(13), iflux, i,matching,isom,MOswitch
      parameter (ncomp = 1)
      double precision Pcm,myresult(15),myerror(15),massth,renscale,result(10),error(10)
      double precision result23(ncomp),error23(ncomp),intresdip(ncomp),errorDip(ncomp)
      double precision gauge,gauge2,vrel,kappa,color_resum1,color_resum8
      double complex stsT2xx_Resum
      character slhafilename*200

      QScale = renscale
      Qscalealphas = Qscale
 
! FPS=1d0 for full polarization sum (lightcone gauge) and FPS=0d0 for ghosts (Feynman gauge)
      FPSk1 = gauge
      FPSk2 = gauge2

** flags assignment
c number of Sfermion1 (squark)
      isf1 = iflag(1)
	
c type of initial state partical (squark)
      itt1 = iflag(2)
	
! number of family of initial state partical (squark)
      igen1 = iflag(3)
	
! number of Sfermion2 (antisquark)
      isf2 = iflag(4)

! type of initial state anti-partical (antisquark)
      itt2 = iflag(5)

! number of family of initial state anti-partical (antisquark)
      igen2 = iflag(6)

! imass = 0- take masses from MicrOmegas, 1-diagonalize here	
      imass = iflag(7)

! iflux = 0 - calculate v.sigma, 1- standard cross-section sigma	
      iflux = iflag(8)

! itree = 0- full one loop, 1- only tree level (quicker)	
      itree = iflag(9)

! include Sommerfeld - 0 = no and 1 = yes 
      isom = iflag(10)

! renormalization scheme
      ischeme = iflag(11)

! legacy
      MOswitch = iflag(12)

! choosesol 
      choosesol = iflag(13)

! ===== Set model parameters =====
! reading the MSSM parameters from micrOmegas and setting the masses & mixings & all shortforms e.g. TB,SA
      call DMNLO_ModelPara(MOswitch,slhafilename)
      call DMNLO_ModelIni(slhafilename)
!	call DMNLO_ModelDigest


! Define the CMS \sqrt{s} and the scale Q2
      sqrtS  = dsqrt(Pcm**2 + MSf2(isf1,itt1,igen1)) + dsqrt(Pcm**2 + MSf2(isf2,itt2,igen2))
      sman = sqrtS**2

! ===== Renormalization Scales for Loop integrals =====
! Divergent UV & IR poles for Loop integrals
      UVdiv = 0d0
      IRdiv = 0d0  ! If LoopTools is used IRdiv has always to be set to zero. 
      
! Renormalization scale \mu & switch between dim.reg and mass regularization for IR divergence - xiIR
! xiIR enables a check with mass regularization if set to 0d0 -> DEFAULT is xiIR = 1d0 !!!! 
      xiIR = 1d0
      muSc = QScale

! EpsPoles = -2,-1,0,1,... steers LoopFunctions
! e.g. EpsPoles = -2/-1 lets all LoopFunctions return the IR divergent coefficient corresponding to double/single pole
! checking UV divergence requires setting EpsPole=1 and UVdiv=1 (loop amplitude should be zero or extremely small)
      EpsPole = 0d0

! ===== Definition of Renormalization Scheme, Scale & Input =====
      call Init_RenScheme(ischeme) 
      if(renfail.eq.0) stop

! flux factor (in mykinematics.h)
      fluxflag = iflux

! Generic couplings (here the type and generation is set !!! )
      call stsT2xx_SetCouplings
      call stsT2xx_SetCounterTerms

! FA kinematics & reduced masses
      call stsT2xx_SetKinematics

! Set lower integration limits for dipole method
      x0IEISg = 2*m1**2/(sman-2*m1**2) ! lower limit on x0 for IE-IS gluon
      x0FEISg = 0.95d0 ! lower limit on x FE-IS gluon


! Set lower integration limits for dipole method
      x0IEISQ = 2*m1**2/(sman-2*m1**2) ! lower limit on x for IE-IS quark
      x0FEISQ = 0.95d0 ! lower limit on x for FE-IS quark
	
*******************************************
! Mass threshold
      massth = sqrtS*(muu3+muu4)

      if (sqrtS.le.massth) then
!
      result = 0d0
      result23 = 0d0
      
      else 

! Compute integrated cross-section
! init LoopTools
      call ltini
      call setmudim(muSc**2)
 
      call stsT2xx_IntCS(result,error,itree)
! exit LoopTools
      call ltexi

      if (itree.eq.0) then
      
      call stsT2xx_Bremsstrahl23(result23, error23)

      call stsT2xx_M2DipolePlusDist(intresdip, errorDip)

      endif
      endif

!
! Copy to output vectors
      myresult = 0d0
      if (itree.eq.0) then

      do i=1,10
            myresult(i) = result(i)
            myerror(i) = error(i)
      enddo

      myresult(11) = sum(result(4:6))+sum(result(8:10))
      myerror(11) = dsqrt(sum(error(4:6)**2)+sum(error(8:10)**2))

! real emission
      myresult(12) = result23(1)
      myerror(12) = error23(1)

! plus distribution
      myresult(13) = intresdip(1)
      myerror(13) = errorDip(1)

      else

      myresult(1) = sum(result(1:3))
      myresult(2) = result(7)

      endif

      if(isom.eq.1) then
      if(itree.eq.0) then
! Sommerfeld enhancement matched to the full NLO calculation and evaluated at the resummation scale (in most cases the Bohr scale)
            matching = 0  !0= NLO-subtraction from NLO-Coulomb, 1= NLO-subtraction from full NLO calculation !The difference is just the scale 1 is 
!                                                 the scale used here, 0 uses the scale from resummation (most likly the bohrscale)
      else 
            matching = 2 ! For Sommerfeld enhancement alone 
      endif
      color_resum1 = -CF  ! attractive
      color_resum8 = 1d0/(2*Nc) ! repulsive

      myresult(14) = stsT2xx_Resum(matching,color_resum1,m1,m1,WSf(isf1,itt1,igen1),WSf(isf2,itt2,igen2),sqrtS,color_resum1)*result(1) 
     &    + stsT2xx_Resum(matching,color_resum8,m1,m1,WSf(isf1,itt1,igen1),WSf(isf2,itt2,igen2),sqrtS,color_resum1)*(result(2)+result(3)+result(7))   
     
! Velocity enhanced part of NLO correction
      vrel = kappa(sman,m1**2,m1**2)/(sman-2*m1**2) 
      myresult(15) = -(Alfas*Pi)/vrel*(color_resum1*result(1) + color_resum8*(result(2)+result(3)+result(7)))
      endif

      end
**************************************************************

      subroutine stsT2xx_IntCS(result,error,mytree)

      implicit none

#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      integer c, neval, fail, ncomp, mytree
      double precision result(10), error(10)

      double precision costhmin, costhmax, relaccuracy, absaccuracy
      double precision kappa

      external stsT2xx_DiffCS, stsT2xx_TreeDiffCS, Patterson
      ncomp = 10
! Integration limits for cos(theta)
      costhmin = -1d0
      costhmax = 1d0

! Set requested accuracy
      absaccuracy = 1d-12
      relaccuracy = 1d-5
    
! Full calculation
	    if (mytree.eq.0) then

	    call Patterson(ncomp, costhmin, costhmax, stsT2xx_DiffCS, 
     &    relaccuracy, absaccuracy, neval, fail, result, error)

! Only tree-level calculation
	    else

	    call Patterson(ncomp, costhmin, costhmax, stsT2xx_TreeDiffCS,
     &    relaccuracy, absaccuracy, neval, fail, result, error)

	    endif

c	if( fail .ne. 0 )
c     &    print *, "Failed to reach the desired accuracy."
	    end

*************************************

      subroutine stsT2xx_DiffCS(result, costh)
      implicit none

#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      double precision costh, intfac, intfacMO, avgfac
      double precision result(10)
      double precision kappa,stsT2gg_M2Dipole,M2Color(3)
      double complex stsT2gg_M2Tree
      double precision stsT2gg_M2TreeSinglet,stsT2gg_M2TreeSOctet,stsT2gg_M2TreeAOctet

      double precision stsT2QQbar_M2Dipole, stsT2QQbar_M2Vertex2qG, stsT2QQbar_M2Vertex2sqG, stsT2QQbar_M2Prop
      double precision stsT2QQbar_M2VertexCounter2qG, stsT2QQbar_M2VertexCounter2sqG, stsT2QQbar_M2CountertermsProp
      double precision stsT2QQbar_M2Box,stsT2QQbar_M2Tree
      double precision stsT2QQbar_M2TreeDipole,Nf

! Mandelstam variables t=(p1-k1)^2 , u=(p2-k1)^2
! Scattering angle theta between p1 and k1

      tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) - 
     &            0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) - 
     &            0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

      uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
      ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

! Integration factor for 2->2 phase-space integration in CMS 
! (the azimuthal integration constant 2pi is inlcuded in intfac)

      intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)

      intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

! Flux factors for different units of cross-section
	if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
	flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
 	
      else
	
! in units of cm^3/sec
! flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
	flux = 1d0/sqrtS**2
	endif

      result = 0d0
! Number of massless quarks
      Nf = 4d0

! ***************************************
! Contributions from stsT2gg
! ***************************************
! Average factor for colours (1/Nc^2) and Bose symmetry factor for identical particles in the final state (1/2)
      avgfac = 1d0/Nc**2/2d0 

! Tree level
!      result(1) = intfac*flux*avgfac*stsT2gg_M2Tree()
! singlet - tree
      result(1) = intfac*flux*avgfac*(stsT2gg_M2TreeSinglet())

! symmetric octet - tree
      result(2) = intfac*flux*avgfac*(stsT2gg_M2TreeSOctet())
! antisymmetric octet - tree
      result(3) = intfac*flux*avgfac*(stsT2gg_M2TreeAOctet())

! 	Dipole contribution
      result(4) = intfac*flux*avgfac*stsT2gg_M2Dipole()

! -------------------------------------   Corrections ----------------------------------------     
!	Vertex contributions
      call stsT2gg_calc_M2Vertex4(M2Color)
	result(5) = intfac*flux*avgfac*sum(M2Color)

      call stsT2gg_calc_M2Vertex3g(M2Color)
      result(5) = result(5) + intfac*flux*avgfac*sum(M2Color)
     
      call stsT2gg_calc_M2Vertex2sg(M2Color)
      result(5) = result(5) + intfac*flux*avgfac*sum(M2Color)

      call stsT2gg_calc_M2VertexGh(M2Color) 
      result(5) = result(5) -intfac*flux*avgfac*(1d0-FPSk1)*(1d0-FPSk2)*sum(M2Color)   

!     Ghost box contributions
      call stsT2gg_calc_M2GhostBoxes(M2Color)
      result(5) = result(5) -intfac*flux*avgfac*(1d0-FPSk1)*(1d0-FPSk2)*sum(M2Color)

!	Propagator contributions
      call stsT2gg_calc_M2Propagator(M2Color)
	result(5) = result(5) + intfac*flux*avgfac*sum(M2Color)

! -------------------------------------   Counterterms ---------------------------------------- 
!	Vertex counterterms
      call stsT2gg_calc_M2Vertex4Counter(M2Color)
	result(6) = intfac*flux*avgfac*sum(M2Color)

      call stsT2gg_calc_M2Vertex2sgCounter(M2Color)
      result(6) = result(6) + intfac*flux*avgfac*sum(M2Color)

      call stsT2gg_calc_M2Vertex3gCounter(M2Color)
      result(6) = result(6) + intfac*flux*avgfac*sum(M2Color)

      call stsT2gg_calc_M2VertexGhCounter(M2Color)
      result(6) = result(6) - intfac*flux*avgfac*(1d0-FPSk1)*(1d0-FPSk2)*sum(M2Color)

!	Propagator counterterms
      call stsT2gg_calc_M2PropagatorCounter(M2Color)
	result(6) = result(6) + intfac*flux*avgfac*sum(M2Color)

      

! ***************************************
! Contributions from stsT2QQbar
! ***************************************
      
! Averaging factor for colour
      avgfac = 1d0/Nc**2
! Tree level
      result(7) = Nf*intfac*flux*avgfac*stsT2QQbar_M2Tree()

! dipoles
      result(8) = intfac*flux*avgfac*stsT2QQbar_M2Dipole()

! Vertex,propagator and box corrections
      result(9) = intfac*flux*avgfac*(stsT2QQbar_M2Vertex2qG() + stsT2QQbar_M2Vertex2sqG() + stsT2QQbar_M2Prop() + stsT2QQbar_M2Box())

! Vertex and propagator counterterms
      result(10) = intfac*flux*avgfac*(stsT2QQbar_M2VertexCounter2qG()+ stsT2QQbar_M2VertexCounter2sqG() + stsT2QQbar_M2CountertermsProp())

      end

      subroutine stsT2xx_TreeDiffCS(result, costh)

      implicit none

#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      double precision costh, intfac, intfacMO, avgColor,avgBose
      double precision result(10)

      double precision kappa,stsT2QQbar_M2Tree,Nf,stsT2gg_M2TreeSinglet,stsT2gg_M2TreeSOctet,stsT2gg_M2TreeAOctet
      double complex stsT2gg_M2Tree

! Mandelstam variables t=(p1-k1)^2 , u=(p2-k1)^2
! Scattering angle theta between p1 and k1

      tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) - 
     &            0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) - 
     &            0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

      uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
      ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

! Integration factor for 2->2 phase-space integration in CMS 
! the azimuthal integration constant 2pi is inlcuded in intfac
      intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)

      ! Flux factor
      if (fluxflag.eq.1) then
! in units of pb
cc	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
      flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else

! in units of cm^3/sec
!     flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
      flux = 1d0/sqrtS**2
      endif

      Nf = 4d0
      avgColor=1d0/Nc**2 
      avgBose = 1d0/2d0

      result = 0d0

! singlet - tree
      result(1) = intfac*flux*avgColor*avgBose*stsT2gg_M2TreeSinglet()
! symmetric octet - tree
      result(2) = intfac*flux*avgColor*avgBose*stsT2gg_M2TreeSOctet()
! antisymmetric octet - tree
      result(3) = intfac*flux*avgColor*avgBose*stsT2gg_M2TreeAOctet()
      result(7) = Nf*intfac*flux*avgColor*stsT2QQbar_M2Tree()

      end

      subroutine stsT2xx_M2DipolePlusDist(intresdip, error)
      implicit none 
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

      
! ***** Variables used for the integration method *****

      integer ndim, ncomp
      parameter (ndim = 2)
      parameter (ncomp = 1)


      integer flags, mineval, maxeval, key, nregions, neval, fail
      integer maxiter, nstart, nincrease

      double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
      integer failcuh
      double precision error(ncomp), prob(ncomp)
      double precision intresdip(ncomp)

      double precision RelAccuracy, AbsAccuracy,kappa
      double precision xi(ndim), fi(ncomp)

      external vegas 
      external stsT2xx_M2DipolePlus

!******************* Integration *********************

! ===== Settings =====

      flags = 0
      MinEval = 1000
      MaxEval = 100000
      MaxIter = 1000
      Nstart = MaxEval/MaxIter
      Nincrease = MaxEval/(2*MaxIter)

! ===== Set desired accuracy =====

      RelAccuracy = 5D-3
      AbsAccuracy = 1D-16

      intresdip = 0d0
      error = 0d0

      if (EpsPole.eq.0d0) then


! Flux factors for different units of cross-section
      if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
      flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else
            
! in units of cm^3/sec
!	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2	 
      flux = 1d0/sqrtS**2
      endif


! ===== VEGAS integration =====

      call vegas(ndim, ncomp, stsT2xx_M2DipolePlus,
     &      dble(RelAccuracy), dble(AbsAccuracy),
     &      flags, MinEval, MaxEval,
     &      Nstart, Nincrease,
     &      Neval, fail, intresdip, error, prob)


      intresdip = intresdip*flux
      error = error*flux

      endif

      end


      subroutine  stsT2xx_M2DipolePlus(ndim, xi, ncomp, fi)
      implicit none     
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      integer ndim, ncomp
      double precision xi(ndim), fi(ncomp) 
      double precision stsT2gg_M2DipolePlusFEIS,stsT2gg_M2DipolePlusIEIS 
      double precision stsT2QQbar_M2DipolePlusFEIS,stsT2QQbar_M2DipolePlusIEIS
      double precision avgColor,avgBose,intfac,Nf

! Number of massless quark flavours
      Nf = 4d0

! ===== Average factor =====
! 1/Nc^2 for color and 1/2 for identical particles in the final state
      avgColor = 1d0/(Nc**2)
      avgBose = 1d0/2d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! stsT2gg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************* FE-IS/IE-FS: Integration for k1/p1 as emitter/spectator pair  ******************* 
      call stsT2xx_SetDipoleKinematicsFEIS(xi,x0FEISg,sman,m1**2,m1**2,0d0,0d0,intfac)

      fi(1) = intfac*avgColor*avgBose*stsT2gg_M2DipolePlusFEIS(1)

!******************* FE-IS/IE-FS: Integration for k1/p2 as emitter/spectator pair  ******************* 

      call stsT2xx_SetDipoleKinematicsFEIS(xi,x0FEISg,sman,m1**2,m1**2,0d0,0d0,intfac)

      fi(1) = fi(1) + intfac*avgColor*avgBose*stsT2gg_M2DipolePlusFEIS(0)

!******************* IE-FS:  p1 as emitter and p2 as spectator  *****************************************
      call stsT2xx_SetDipoleKinematicsIEIS(xi,x0IEISg,sman,m1**2,m1**2,0d0,0d0,intfac) 

      fi(1) = fi(1) + intfac*avgColor*avgBose*stsT2gg_M2DipolePlusIEIS()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! stsT2QQbar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************* FE-IS/IE-FS: Integration for k1/p1 as emitter/spectator pair  ******************* 
      call stsT2xx_SetDipoleKinematicsFEIS(xi,x0FEISQ,sman,m1**2,m1**2,0d0,0d0,intfac)

      fi(1) = fi(1) + intfac*avgColor*Nf*stsT2QQbar_M2DipolePlusFEIS(1)

!******************* FE-IS/IE-FS: Integration for k1/p2 as emitter/spectator pair  ******************* 

      call stsT2xx_SetDipoleKinematicsFEIS(xi,x0FEISQ,sman,m1**2,m1**2,0d0,0d0,intfac)

      fi(1) = fi(1) + intfac*avgColor*Nf*stsT2QQbar_M2DipolePlusFEIS(0)

!******************* IE-FS:  p1 as emitter and p2 as spectator  *****************************************
      call stsT2xx_SetDipoleKinematicsIEIS(xi,x0IEISQ,sman,m1**2,m1**2,0d0,0d0,intfac) 

      fi(1) = fi(1) + intfac*avgColor*Nf*stsT2QQbar_M2DipolePlusIEIS()

      end




      subroutine stsT2xx_Bremsstrahl23(intres, error)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      integer ndim, ncomp, iter
      parameter (ndim = 4)
      parameter (ncomp = 1)
      integer flags, mineval, maxeval, key, nregions, neval, fail
      integer maxiter, nstart, nincrease, nnew

      double precision intfac23
      double precision RELACCURACY, ABSACCURACY
      double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
      integer failcuh
      double precision error(ncomp), prob(ncomp)
      double precision intres(ncomp)
      double precision kappa


      external vegas
      external cuhre
      external suave
      external mint
      external stsT2xx_Integrand4
                 
                  
      flags = 0
      MinEval = 1000
      MaxEval = 100000
      MaxIter = 1000
      Nstart = MaxEval/MaxIter
      Nincrease = MaxEval/(2*MaxIter)
      
      RELACCURACY = 5D-3
      ABSACCURACY = 1D-16
                  

      call vegas(ndim, ncomp, stsT2xx_Integrand4,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      flags, MINEVAL, MAXEVAL,
     &      nstart, nincrease,
     &      neval, fail, intres, error, prob)
      
                  
c integration factor already includes the trivial integration over phi
      
      intfac23 = sqrtS**2/(32d0*(2*pi)**4)
      
c* flux 
      
      if (fluxflag.eq.1) then		!Mo: neu von mir hinzugefuegt, jetzt kann auch der Prozess 2->3 ohne den Prozess 2->2 betrachtet werden, da hier 
c* in units of pb                                                                       !    jetzt des Flußfaktor unabhaengig nochmal bestimmt wird.
cc	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
c* in units of GeV^-2
      flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else
      
c* in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
c* in units of GeV^-2
      flux = 1d0/sqrtS**2
      endif

                             
      intres = flux*intfac23*intres
      error = flux*intfac23*error	
      
      end

      subroutine stsT2xx_Integrand4(ndim, x, ncomp, f)
      implicit none

#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

      
      integer ndim, ncomp
      double precision x(*), f(*)

      double precision cosxi,avgColor,avgBose,Nf
      double precision stsT2xx_M2RealQ, stsT2xx_M2DipoleRealQ
!      double complex stsT2gg_M2Real3gPhysical
      double precision stsT2xx_M2Real3g,stsT2xx_M2RealGhost,stsT2xx_M2DipoleReal3g
      integer io_error

! ===== Integration boundaries and Jacobi determinant =====

      call stsT2xx_PhSpace(x,muu3,muu4,muu5,muu1,muu2,cosxi)      
      
      f(1) = 0d0

! Number of massless quark flavours
      Nf = 4d0

! ===== Single contributions to matrix element =====

      if(dabs(cosxi).lt.0.999999999999999d0) then              !cut out numerical errors leading to a |cosxi|>1 oder cosxi=NAN (tests x1 and x2)

! 1/Nc^2 for color and 1/6 for identical particles in the final state
      avgColor = 1d0/(Nc**2)
      avgBose = 1d0/6d0

! Gluon contribution 
      f(1) = f(1) + jacobi*avgColor*avgBose*(stsT2xx_M2Real3g() -  stsT2xx_M2RealGhost() - stsT2xx_M2DipoleReal3g())

! Quark contribution 
      f(1) = f(1) + Nf*jacobi*avgColor*(stsT2xx_M2RealQ() - stsT2xx_M2DipoleRealQ())

      endif
      
      end	
                    
               
! ================= General Phase-space element======================
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine stsT2xx_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb,cosxi)
      implicit none
              
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      double precision x(*)
      double precision mmr1,mmr2,mmr3,mmra,mmrb
      double precision eta, costh
      double precision cosxi,sinxi,sinth
      double precision xx3max,xx3min,xx1max,xx1min
      double precision taubar,sigmabar,muplus,muminus

      xx3min = 2d0*mmr3/sqrtS   
      xx3max = 1d0 - (mmr1 + mmr2)**2 + mmr3**2

      xx3 = (xx3max - xx3min)*x(2) + xx3min
     
      taubar = 1d0 - xx3 + mmr3**2
      sigmabar = 2d0 - xx3
      muplus = mmr1 + mmr2
      muminus = mmr1 - mmr2

      xx1min = 1d0/(2d0*taubar)*(sigmabar*(taubar+muplus*muminus) - dsqrt(xx3**2-4d0*mmr3**2)*dsqrt((taubar-muplus**2)*(taubar-muminus**2))) 
      xx1max = 1d0/(2d0*taubar)*(sigmabar*(taubar+muplus*muminus) + dsqrt(xx3**2-4d0*mmr3**2)*dsqrt((taubar-muplus**2)*(taubar-muminus**2)))
      xx1 = (xx1max - xx1min)*x(1) + xx1min

      xx2 = 2d0 - xx1 - xx3

      eta = 2*pi*x(3)
      costh = 2*x(4)-1d0

      jacobi = (xx1max-xx1min)*(xx3max-xx3min)*2d0*pi*2d0
      
      cosxi = ((2 - xx1 - xx3)**2 + 4*mmr1**2 + 4*mmr3**2 - 4*mmr2**2 - xx1**2 - xx3**2)/(2d0*dsqrt(xx1**2 - 4*mmr1**2)*dsqrt(xx3**2 - 4*mmr3**2))

      sinxi = dsqrt(1-cosxi**2)
      sinth = dsqrt(1-costh**2)

 ! *******************************Scalar products ******************************                                                                                                                                   
      p1p2 = (sqrtS**2/2d0)*(1d0-mmra**2-mmrb**2)
      p1k1 = sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
      p1k3 = sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)                                                                                       
      p1k2 = (sqrtS**2/2d0)*(1d0-mmrb**2+mmra**2)-p1k1-p1k3				!Momentumconservation: p1k2=p1p1+p1p2-p1k1-p1k3
       
      k1k2 = (1d0/2d0)*sqrtS**2*(1d0 - xx3 - mmr1**2 - mmr2**2 + mmr3**2)
      k1k3 = (1d0/2d0)*sqrtS**2*(1d0 - xx2 - mmr1**2 + mmr2**2 - mmr3**2)
      k2k3 = (1d0/2d0)*sqrtS**2*(1d0 - xx1 + mmr1**2 - mmr2**2 - mmr3**2)

      p2k1 = -p1k1+sqrtS**2*mmr1**2+k1k2+k1k3
      p2k3 = -p1k3+sqrtS**2*mmr3**2+k1k3+k2k3
      p2k2 = -p1k2+sqrtS**2*mmr2**2+k1k2+k2k3

! These scalar products are necessary for the physical gauge as they enter through the polarization sum.
! When using ghosts, these scalar products can be commented out.
! The "bar" corresponds to a parity transformed four-vector. 
!      p1Barp2 = ((-mmra**4 + mmrb**2 - mmrb**4 + mmra**2*(1 + 2*mmrb**2))*sqrtS**2)/2d0
!      p1Bark1 = ((1 + mmra**2 - mmrb**2)*sqrtS**2*xx1)/4. + (dsqrt(mmra**4 + (-1 + mmrb**2)**2 - 2*mmra**2*(1 + mmrb**2))*sqrtS**2*dsqrt(-4*mmr1**2 + xx1**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))/4d0
!      p1Bark3 =(sqrtS**2*((1 + mmra**2 - mmrb**2)*xx3 + costh*dsqrt(mmra**4 + (-1 + mmrb**2)**2 - 2*mmra**2*(1 + mmrb**2))*dsqrt(-4*mmr3**2 + xx3**2)))/4d0                                                           
!      p1Bark2 = (sqrtS**2/2d0)*(1d0-mmrb**2+mmra**2)-p1Bark1-p1Bark3				
!
!      k1Bark2 = (sqrtS**2*(-xx1**2 + xx1*(2 + xx2 - xx3) + 2*(-1 + mmr1**2 + mmr2**2 - mmr3**2 + xx3)))/4d0
!      k1Bark3 = (sqrtS**2*(1 + mmr1**2 - mmr2**2 + mmr3**2 + xx1*(-1 + xx3) - xx3))/2d0
!      k2Bark3 = (sqrtS**2*(2*(-1 - mmr1**2 + mmr2**2 + mmr3**2) - xx1*(-2 + xx3) + (2 + xx2)*xx3 - xx3**2))/4d0
!
!      
!      k1Bark1=(sqrtS**2*(-2*mmr1**2 + xx1**2))/2d0
!      k2Bark2=(sqrtS**2*(4 - 4*mmr2**2 + xx1**2 + xx2**2 + 2*xx1*(-2 + xx3) - 4*xx3 + xx3**2))/4d0
!      k3Bark3=(sqrtS**2*(-2*mmr3**2 + xx3**2))/2d0
!      p1Barp1=((mmra**4 - 2*mmra**2*mmrb**2 + (-1 + mmrb**2)**2)*sqrtS**2)/2d0
!      p2Barp2=((1 + mmra**4 + mmrb**4 - 2*mmra**2*(1 + mmrb**2))*sqrtS**2)/2d0
!
!      p2Bark1 = -p1Bark1+k1Bark1+k1Bark2+k1Bark3
!      p2Bark3 = -p1Bark3+k1Bark3+k2Bark3+k3Bark3
!      p2Bark2 = -p1Bark2+k1Bark2+k2Bark2+k2Bark3

      end


#include "common/stsT2xx_Couplings.F"
#include "common/stsT2xx_Kinematics.F"
#include "common/stsT2xx_SetCounterterms.F"
#include "common/stsT2xx_Loops.F"
#include "common/stsT2xx_DipoleKinematics.F"

#include "common/stsT2xx_Resum.F"
#include "common/stsT2xx_LoopTools.F"


#include "stsT2gg/stsT2gg_M2Tree.F"
#include "stsT2gg/stsT2gg_M2Vertex4.F"
#include "stsT2gg/stsT2gg_M2Vertex3g.F"
#include "stsT2gg/stsT2gg_M2Vertex2sg.F"
#include "stsT2gg/stsT2gg_M2Counterterms.F"
#include "stsT2gg/stsT2gg_M2Propagator.F"
#include "stsT2gg/stsT2gg_M2Ghost.F"
#include "stsT2gg/stsT2gg_M2Dipole.F"
#include "stsT2gg/stsT2gg_M2DipolePlus.F"


#include "stsT2QQbar/stsT2QQbar_M2Box.F"
#include "stsT2QQbar/stsT2QQbar_M2Counterterms2qG.F"
#include "stsT2QQbar/stsT2QQbar_M2Counterterms2sqG.F"
#include "stsT2QQbar/stsT2QQbar_M2CountertermsProp.F"
#include "stsT2QQbar/stsT2QQbar_M2Dipole.F"
#include "stsT2QQbar/stsT2QQbar_M2DipolePlus.F"
#include "stsT2QQbar/stsT2QQbar_M2Prop.F"
#include "stsT2QQbar/stsT2QQbar_M2Tree.F"
#include "stsT2QQbar/stsT2QQbar_M2Vertex2qG.F"
#include "stsT2QQbar/stsT2QQbar_M2Vertex2sqG.F"


#include "stsT2xx_M2Real3g.F"
#include "stsT2xx_M2RealGhost.F"
#include "stsT2xx_M2RealQ.F"
#include "stsT2xx_M2DipoleReal.F"
! It is advised to use the physical gauge only for testing purposes due to much longer compilation and runtimes.
!#include "stsT2xx_M2Real3gPhysical.F"


