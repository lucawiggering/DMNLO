      subroutine staustau2QQ(renscale,Pcm,iflag,myresult,slhafilename)
      implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "../util/DMNLO_RenScheme.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_GenCouplings.h" 

      integer ndim, ncomp, iflag(17), inptype, inpgen, iflux,i,MOswitch
      parameter (ncomp = 1)

      double precision Pcm,myresult(9), sinth, costh, avgfac, massth,renscale
      double precision result(9), error(9)
      double precision result23(ncomp),error23(ncomp),result23num(ncomp),error23num(ncomp)
      double complex staustau2QQ_M2Brems,staustau2QQ_M2enhancement
      character slhafilename*200
	
      QScale = renscale
      Qscalealphas = Qscale

** flags assignment
c sfermion index of incoming particle
      isf1 = iflag(1)
	
c type of incoming particle
      itt1 = iflag(2)
	
! generation of incoming particle
      igen1 = iflag(3)
	
c sfermion index of incoming anti-particle
      isf2 = iflag(4)
		
! type of incoming anti-particle
      itt2 = iflag(5)

! generation of incoming anti-particle
      igen2 = iflag(6)
	
! type of final state particle
      ftt1 = iflag(7)

! generation of final state particle
      fgen1 = iflag(8)
	
! type of final state anti-particle
      ftt2 = iflag(9)

! generation of final state anti-particle
      fgen2 = iflag(10)

! imass = 0- take masses from MicrOmegas, 1-diagonalize here	
      imass = iflag(11)

! iflux = 0 - calculate v.sigma, 1- standard cross-section sigma	
      iflux = iflag(12)

! itree = 0- full one loop, 1- only tree level (quicker)	
      itree = iflag(13)

! ischeme - choose renormalization scheme
      ischeme = iflag(14)

! isomm = 1 - add Sommerfled enhancement, 0 don't
      isomm = iflag(15) ! in Model.h

! legacy
      MOswitch = iflag(16)

! choosesol
      choosesol = iflag(17)

! ===== Set model parameters =====
! reading the MSSM parameters from micrOmegas and setting the masses & mixings & all shortforms e.g. TB,SA

      call DMNLO_ModelPara(MOswitch,slhafilename)
      call DMNLO_ModelIni(slhafilename)
!	call DMNLO_ModelDigest

      sqrtS  = dsqrt(Pcm**2 + MSf2(isf1,itt1,igen1)) + dsqrt(Pcm**2 + MSf2(isf2,itt2,igen2))

! ===== Renormalization Scales for Loop integrals =====
! Divergent UV & IR poles for Loop integrals
      UVdiv = 0d0
      IRdiv=0d0
	
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
      call staustau2QQ_SetCouplings

! FA kinematics & reduced masses
      call staustau2QQ_SetKinematics

*******************************************

! Mass threshold
      massth = sqrtS*(muu3+muu4)


      if (sqrtS.le.massth) then
      do i=1,9
      result(i) = 0d0
      enddo
      result23(1) = 0d0
      else 

! Compute integrated cross-section

     	call staustau2QQ_IntCS(result,error,itree)

! Copy to output vectors
         
      if (itree.eq.0) then

            call staustau2QQ_Bremsstrahl23(result23num,error23num)

            do i=1,6
            myresult(i) = result(i)
            enddo
            myresult(9) = result(1)+result(2)+result(3)+result(5)+result(6)+result23num(1)
      else

            do i=1,9
            myresult(i) = result(i)
            enddo

      endif
      endif
	end


      subroutine staustau2QQ_IntCS(result,error,mytree)
      implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"

      double precision result(9), error(9)
      double precision costhmin, costhmax, relaccuracy, absaccuracy
      double precision kappa
      integer c, neval, fail, ncomp, mytree
      external staustau2QQ_TreeDiffCS, staustau2QQ_DiffCS, Patterson

      ncomp = 7

! Integration limits for cos(theta)
      costhmin = -1d0
      costhmax = 1d0

! Set requested accuracy
      absaccuracy = 1d-12
      relaccuracy = 1d-5

! Full calculation
      if (mytree.eq.0) then

      call Patterson(ncomp, costhmin, costhmax, staustau2QQ_DiffCS, 
     &    relaccuracy, absaccuracy, neval, fail, result, error)

! Only tree-level calculation
      else

      call Patterson(ncomp, costhmin, costhmax, staustau2QQ_TreeDiffCS,
     &    relaccuracy, absaccuracy, neval, fail, result, error)

      endif

	end

 
      subroutine staustau2QQ_DiffCS(result, costh)

      implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac, testmin, testmax
      double precision result(9)
      double precision kappa,staustau2QQ_M2Dipole
      double complex staustau2QQ_M2Tree, staustau2QQ_M2Vertex, staustau2QQ_M2Counter,staustau2QQ_M2enhancement,staustau2QQ_M2TreeH, staustau2QQ_M2TreeV
	
! Mandelstam variables t=(p1-k1)^2 , u=(p1-k2)^2
      tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

      uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
      ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred


! Integration factor for 2->2 phase-space integration in CMS 
! (the azimuthal integration constant 2pi is inlcuded in intfac)

      intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)

! Flux factors for different units of cross-section
      if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
            flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else
	
! in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif
!!distinguish between identical and non identical particels
!	if(indk.ge.indl) then 
!		if(indk.lt.5)then
!! Average factor+Symmetry factor from equal final particels
!	avgfac = 1d0/4d0*1d0/2d0
!		end if
!	else
! Average factor
      avgfac = 1d0
!	end if
! ***************************************
! contributions to the matrix element 
! ***************************************
      result = 0d0
      result(1) = intfac*flux*avgfac*staustau2QQ_M2Tree()
      result(2) = intfac*flux*avgfac*staustau2QQ_M2Vertex()
      result(3) = intfac*flux*avgfac*staustau2QQ_M2Counter()
      if(isomm.eq.1) then	  
            result(5) = staustau2QQ_M2enhancement()*result(1)    
      endif
      result(6) = intfac*flux*avgfac*staustau2QQ_M2Dipole()
	end

*************************************

      subroutine staustau2QQ_TreeDiffCS(result, costh)

      implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
      double precision costh, intfac, avgfac,result(9),kappa
      double complex staustau2QQ_M2Tree,staustau2QQ_M2enhancement

! Mandelstam variables t=(p2-k1)**2 , u=(p1-k1)**2
! Scattering angle theta between p1 and k1

      tman = sqrtS**2*(muu2**2 + muu3**2 - 0.5d0*(1d0-muu1**2+muu2**2)*(1d0-muu4**2+muu3**2) - 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu2**2 + muu3**2 - 0.5d0*(1d0-muu1**2+muu2**2)*(1d0-muu4**2+muu3**2) - 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

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
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif
      avgfac = 1d0
	
! no flux and averaging factor!
      if(isomm.ne.1) then
            result(1) = intfac*flux*avgfac*staustau2QQ_M2Tree()
      else 
            result(1) = (staustau2QQ_M2enhancement() + 1d0)*intfac*flux*avgfac*staustau2QQ_M2Tree() !!factorised LO amplitude with enhancement (Green's function) added to LO
      endif

	  
	end

      subroutine staustau2QQ_Bremsstrahl23(intres, error)

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"

      integer ndim, ncomp
      parameter (ndim = 4)
      parameter (ncomp = 1)
      integer flags, mineval, maxeval, key, nregions, neval, fail
      integer maxiter, nstart, nincrease

      double precision intfac23, avgfac, x1min, x1max, x2min, x2max
      double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
      integer failcuh
      double precision error(ncomp), prob(ncomp)
      double precision intres(ncomp)
      double precision kappa

      external integrand4
      external staustau2QQ_integrand4
      external vegas
      external cuhre


	flags = 0
	mineval = 1000
	maxeval = 50000
	maxiter = 100
	nstart = maxeval/MAXITER
	nincrease = maxeval/(2*MAXITER)

	RELACCURACY = 5D-3
	ABSACCURACY = 1D-16

        	  call vegas(ndim, ncomp, staustau2QQ_integrand4,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      flags, MINEVAL, MAXEVAL,
     &      nstart, nincrease,
     &      neval, fail, intres, error, prob)


 
! integration factor already includes the trivial integration over phi

	intfac23 = sqrtS**2/(32d0*(2*pi)**4)

!* flux 

      if (fluxflag.eq.1) then																!Mo: neu von mir hinzugefuegt, jetzt kann auch der Prozess 2->3 ohne den Prozess 2->2 betrachtet werden, da hier 
!* in units of pb                                                                       !    jetzt des Flußfaktor unabhaengig nochmal bestimmt wird.
!c	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
!* in units of GeV^-2
            flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else

!* in units of cm^3/sec
!	 flux = 2.9979d-26*hbar_c2/sqrtS**2
!* in units of GeV^-2
            flux = 1d0/sqrtS**2
	endif

      avgfac = 1d0/4d0                                                                

      intres(1) = avgfac*flux*intfac23*intres(1)
      error(1) = avgfac*flux*intfac23*error(1)

	end
	
!====================================================================
      subroutine staustau2QQ_Integrand4(ndim, x, ncomp, f)

      implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
      integer ndim, ncomp
      double precision x(*), f(*)

      double precision eta,costh,x1max,x1min,x2max,x2min
      double precision cosxi, sinxi, sinth
      double precision staustau2QQ_M2Radiation

      external staustau2QQ_PhSpace

! ===== Integration boundaries and Jacobi determinant =====

! dipole configuration -> gluon integrated out

      call staustau2QQ_PhSpace(x,muu3,muu5,muu4,muu1,muu2,cosxi)            !exchange masses as needed for correct 2->3 phase-space integration; integration over Quark1 and Quark2 phase space (for quark-quark-Gluon-finalstate)

! ===== Scalar products - dipole configuration =====

      x1 = xx1
      x2 = xx3    

	
! ===== Single contributions to matrix element =====

      if(dabs(cosxi).lt.0.999999999999999d0) then              !cut out numerical errors leading to a |cosxi|>1 oder cosxi=NAN (tests x1 and x2)

      f(1) = 0d0					    !Die Wirkungsquerschnitte in Radiation-Datei enthalten integriert Anteile
                                                      !=> Ausgabe des differentiellen WQ's also nicht mehr möglich !!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      f(1) = jacobi*staustau2QQ_M2Radiation()  

      else

      f(1) = 0d0

      endif


	end

! ================= General Phase-space element======================
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine staustau2QQ_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb,cosxi)

      implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"

      double precision x(*)
      double precision mmr1,mmr2,mmr3,mmra,mmrb
      double precision eta, costh
      double precision cosxi,sinxi,sinth
      double precision xx3max,xx3min,xx1max,xx1min
      double precision taubar,sigmabar,muplus,muminus


! ***** Integration boundaries and Jacobi determinant (from HitchHiker) *****

      xx3min = 2*mmr3              !k3->k2

                  
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

! *******************************Scalar products *******************************


!Diple subtraction

       p1p2 = (sqrtS**2/2d0)*(1d0-mmra**2-mmrb**2)
       p1k1 = sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
       p1k3 = sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)
       p1k2 = (sqrtS**2/2d0)*(1d0-mmrb**2+mmra**2)-p1k1-p1k3             !Momentumconservation: p1k2=p1p1+p1p2-p1k1-p1k3
 
       k1k2 = (1d0/2d0)*sqrtS**2*(1d0 - xx3 - mmr1**2 - mmr2**2 + mmr3**2)
       k1k3 = (1d0/2d0)*sqrtS**2*(1d0 - xx2 - mmr1**2 + mmr2**2 - mmr3**2)
       k2k3 = (1d0/2d0)*sqrtS**2*(1d0 - xx1 + mmr1**2 - mmr2**2 - mmr3**2)

       p2k1 = -p1k1+sqrtS**2*mmr1**2+k1k2+k1k3
       p2k3 = -p1k3+sqrtS**2*mmr3**2+k1k3+k2k3
       p2k2 = -p1k2+sqrtS**2*mmr2**2+k1k2+k2k3         

	end

#include "staustau2QQ_Couplings.F"
#include "staustau2QQ_Kinematics.F"
#include "staustau2QQ_M2Tree.F"
#include "staustau2QQ_M2TreeH.F"
#include "staustau2QQ_M2TreeV.F"
#include "staustau2QQ_M2Counter.F"
#include "staustau2QQ_M2Vertex.F"
#include "staustau2QQ_M2Brems.F"
#include "staustau2QQ_Enhancement.F"
#include "staustau2QQ_M2Dipole.F"
#include "staustau2QQ_M2Radiation.F"
