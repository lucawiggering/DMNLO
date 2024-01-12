
      subroutine stst2qq(renscale,Pcm,iflag,myresult,myerror,slhafilename)
      implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_GenCouplings.h" 
#include "stst2QQ_GenCounterterms.h"

      integer ndim, ncomp, iflag(18),iflux, i,iQCDEW,MOswitch,matching
      parameter (ncomp = 1)
!
      double precision Pcm,myresult(13),massth, myerror(13),renscale
      double precision result(13), error(13)
      double precision result23(ncomp),error23(ncomp)
      double precision intresdip(ncomp), errordip(ncomp)
      double precision result_Coulomb,color_resum3,color_resum6,stst2QQ_Sommerfeld
      character slhafilename*200
      double complex DMNLO_Resum

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

! iQCDEW = 0 - all QCD^2 + QCD-EW + EW^2, 
!          1 - only QCD^2 , 2 - only QCD-EW interference, 3 - only EW^2
!         12 - QCD^2 + QCD-EW, 13 - QCD^2 + EW^2, 23 - QCD-EW + EW^2 
      iQCDEW = iflag(14)

! isomm = 1 - add Sommerfled enhancement, 0 don't
      isomm = iflag(15) ! in Model.h

! ischeme - choose renormalization scheme
      ischeme = iflag(16)

! legacy
      MOswitch = iflag(17)

! choosesol
      choosesol = iflag(18) !Set which solution should be taken for the Mixmatrices USf - only implemented in Renscheme 1 (0:DMNLO chooses solution, 1: take solution1, 2: take solution 2)
      
      Res_Style = 1
      NNLO_Res = 1

! iQCD flag switches between different orders at tree-level
      if (iQCDEW.eq.0) then
            iQCD(1) = 1d0
            iQCD(2) = 1d0
            iQCD(3) = 1d0
      else if (iQCDEW.eq.1) then
            iQCD(1) = 1d0
            iQCD(2) = 0d0
            iQCD(3) = 0d0
      else if (iQCDEW.eq.2) then
            iQCD(1) = 0d0
            iQCD(2) = 1d0
            iQCD(3) = 0d0
      else if (iQCDEW.eq.3) then
            iQCD(1) = 0d0
            iQCD(2) = 0d0
            iQCD(3) = 1d0
            else if (iQCDEW.eq.12) then
                  iQCD(1) = 1d0
                  iQCD(2) = 1d0
                  iQCD(3) = 0d0
            else if (iQCDEW.eq.13) then
                  iQCD(1) = 1d0
                  iQCD(2) = 0d0
                  iQCD(3) = 1d0
            else if (iQCDEW.eq.23) then
                  iQCD(1) = 0d0
                  iQCD(2) = 1d0
                  iQCD(3) = 1d0
      endif


! ===== Set model parameters =====
! reading the MSSM parameters from micrOmegas and setting the masses & mixings & all shortforms e.g. TB,SA

      call DMNLO_ModelPara(MOswitch,slhafilename)
      call DMNLO_ModelIni(slhafilename)
!	call DMNLO_ModelDigest

! ===== Renormalization Scales for Loop integrals =====
! Divergent UV & IR poles for Loop integrals
      UVdiv = 0d0
      IRdiv = 0d0
    
! Renormalization scale \mu & switch between dim.reg and mass regularization for IR divergence - xiIR
! xiIR enables a check with mass regularization if set to 0d0 -> DEFAULT is xiIR = 1d0 !!!! 
      xiIR = 1d0
      muSc = QScale

! EpsPoles = -2,-1,0,1,... steers LoopFunctions
! e.g. EpsPoles = -2/-1 lets all LoopFunctions return the IR diveregent coefficient corresponding to double/single pole
! checking UV divergence requires setting EpsPole=1 and UVdiv=1 (loop amplitude should be zero or extremely small)
      EpsPole = 0d0

! ! ===== Definition of Renormalization Scheme, Scale & Input =====
	call Init_RenScheme(ischeme)			!Pass to 1 again -- 99 for comparison with MO
      if(renfail.eq.0) stop
      
! Computation of CM energy	  
      sqrtS  = dsqrt(Pcm**2 + MSf2(isf1,itt1,igen1)) + dsqrt(Pcm**2 + MSf2(isf2,itt2,igen2))

! flux factor (in mykinematics.h)
      fluxflag = iflux

! Generic couplings (here the type and generation is set !!! )
      call stst2qq_SetCouplings

! FA kinematics & reduced masses
      call stst2qq_SetKinematics

! self-energies & counterterms
      call stst2qq_SetCounterTerms

	
*******************************************

! Mass threshold
      massth = sqrtS*(muu3+muu4)

! Set integration limits for dipole method
      x0IEIS = 2*m1*m2/(sqrtS**2-m1**2-m2**2) ! lower limit for IE-IS
      x0FEIS = 0.95d0 ! lower limit for FE-IS

      if (sqrtS.le.massth) then
      do i=1,9
      result(i) = 0d0
      enddo
      result23(1) = 0d0

      else 

!TODO: check for (isf1.eq.isf2).and.(ftt1.eq.ftt2).and.((isf1.neq.ftt1).or.(isf2.neq.ftt2))  
      call stst2qq_IntCS(result,error,itree)
      
!------------------ For Resummation------------------------------------------------------------------------------------

!setting colorbasis factor - please think about attractive and repulsive when considering the minus sign convention. In DMNLO_Resum is a minus sign assigned to the colorbasis following the convention of the Christoph Borschensky 
      color_resum3 = -(Nc+1d0)/(2*Nc)       !Colorfactor for resummation of triplett
      color_resum6 = (Nc-1d0)/(2*Nc)       !Colorfactor for resummation of sextet

! matching: 0= NLO-subtraction from NLO-Coulomb, 1= NLO-subtraction from full NLO calculation !The difference is just the scale 1 is the scale used here in the full NLO calculation, 0 uses the scale from resummation (most likly the bohrscale), 2= no subtraction of the NLO contribution at all --> Sommerfeld enhancement alone

      matching = 2
      if(itree.eq.0) then
      matching = 0
      result23(1) = 0d0
      call stst2QQ_Bremsstrahl23(result23, error23) 
      result(9) = result23(1)

      intresdip = 0d0
      call stst2QQ_M2DipolePlusDist(intresdip, errordip)
       
      
      endif

! Copy to output vectors
      do i=1,13
      myresult(i) = result(i)
      enddo

      if(isomm.eq.1) then
            myresult(1) =  DMNLO_Resum(matching,0,color_resum3,m1,m2,WSf(isf1,itt1,igen1),WSf(isf2,itt2,igen2),sqrtS)*result(6) 
     &         +  DMNLO_Resum(matching,0,color_resum6,m1,m2,WSf(isf1,itt1,igen1),WSf(isf2,itt2,igen2),sqrtS)*result(7)     
      endif

      myerror(1) = error23(1) ! Bremsstrahl error
      myerror(2) = errordip(1) ! plus-dist-error
      myerror(3) = dsqrt(error(2)**2 + error(3)**2 + error(4)**2 + error(12)**2)  !virtual error (+virtual dipole part)

      myresult(9) = myresult(1) + intresdip(1) + result23(1) + result(2)+ result(3) + result(4) + result(10) + result(12)
! check if NLO cross section is negative, if it is set cross section to 0 
      if(myresult(9).lt.0d0) then
            myresult(9) = 0d0
      endif
      
      endif

      end

**************************************************************

      subroutine stst2qq_IntCS(result,error,mytree)

      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"

      integer ncomp
      parameter (ncomp = 13) 
      double precision result(ncomp), error(ncomp)

      double precision costhmin, costhmax, relaccuracy, absaccuracy
      double precision kappa
      integer c, neval, fail, mytree
      external stst2qq_TreeDiffCS,stst2qq_TreesplitDiffCS, stst2qq_DiffCS, Patterson


! Integration limits for cos(theta)
      costhmin = -1d0
      costhmax = 1d0

! Set requested accuracy
      absaccuracy = 1d-12
      relaccuracy = 1d-5
	 
! 	print*,'mytree',mytree
! Full calculation
      if (mytree.eq.0) then

	  call Patterson(ncomp, costhmin, costhmax, stst2qq_DiffCS, 
     &    relaccuracy, absaccuracy, neval, fail, result, error)
! Only tree-level splittet in seperate contribution calculation     
      elseif (mytree.eq.2) then

	  call Patterson(ncomp, costhmin, costhmax, stst2qq_TreesplitDiffCS, 
     &    relaccuracy, absaccuracy, neval, fail, result, error)
         

! Only tree-level calculation
      else

	  call Patterson(ncomp, costhmin, costhmax, stst2qq_TreeDiffCS,
     &    relaccuracy, absaccuracy, neval, fail, result, error)

	endif

	end

*************************************

      subroutine stst2qq_DiffCS(result, costh)
      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac,result(13),kappa
      double complex stst2qq_M2Tree,stst2QQ_M2Counterterm,stst2QQ_M2CountertermsProp,ststQQ_M2VertexT,stst2QQ_M2Prop
      double complex stst2QQ_M2Boxitjs,stst2QQ_M2Boxjtis,stst2QQ_M2Boxsextet,stst2QQ_M2Boxtriplett
      double complex stst2QQ_M2Boxcancel3,stst2QQ_M2Boxcancel6
      double precision stst2QQ_SoftGluon,stst2QQ_M2Dipole
      double complex stst2QQ_M2realemissionTGLTGL,stst2QQ_M2realemissionUGLUGL,stst2QQ_M2realemissionTGLUGL
      double complex stst2QQ_M2realemissionTNeuUNeu,stst2QQ_M2realemissionTNeuTNeu,stst2QQ_M2realemissionUNeuUNeu
      double complex stst2QQ_M2realemissionTGlUNeu,stst2QQ_M2realemissionUGlTNeu,stst2QQ_M2realemissionTGlTNeu,stst2QQ_M2realemissionUGlUNeu


! Mandelstam variables t=(p1-k1)^2 , u=(p1-k2)^2
      tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
      tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0+muu3**2-muu4**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

      uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
      ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

      vrel = dsqrt(((sqrtS**2-m1**2-m2**2)**2-4*m1**2*m2**2))/(sqrtS**2-(m1**2+m2**2))

! Integration factor for 2->2 phase-space integration in CMS 
! (the azimuthal integration constant 2pi is inlcuded in intfac)
      intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)

! Flux factors for different units of cross-section
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

! Average factor
      if ((ftt1.eq.ftt2).and.(fgen1.eq.fgen2)) then
            avgfac = 1d0/(Nc**2)*1d0/2d0
      else 
            avgfac = 1d0/(Nc**2)
      end if

! ***************************************
! contributions to the matrix element 
! ***************************************
      result = 0d0
      result(1) = intfac*flux*avgfac*stst2QQ_M2Tree()

      result(2) = intfac*flux*avgfac*(ststQQ_M2VertexT()+stst2QQ_M2Prop())
      result(3) = intfac*flux*avgfac*(stst2QQ_M2Counterterm()+stst2QQ_M2CountertermsProp())
      result(4) = intfac*flux*avgfac*(stst2QQ_M2Boxsextet()+stst2QQ_M2Boxtriplett())
      if(isomm.eq.1) then
            result(10) = intfac*flux*avgfac*(-stst2QQ_M2Boxcancel3()-stst2QQ_M2Boxcancel6())
      endif

!      result(5) = intfac*flux*avgfac*stst2QQ_SoftGluon()
      result(12) = intfac*flux*avgfac*stst2QQ_M2Dipole()
!      write(*,*) "IR check",sum(result(2:4))/result(12)
! ! !  ! ! ! ! ! For Resummation- singlet: 
      result(6) = intfac*flux*avgfac*stst2QQ_M2Tree3
! ! ! ! ! For Resummation- octet: 
      result(7) = intfac*flux*avgfac*stst2QQ_M2Tree6

      end
 	  

*************************************

      subroutine stst2qq_TreeDiffCS(result, costh)

      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac
      double precision result(13)


      double precision kappa, stst2qq_M2TreeM0, stst2qq_M2TreePerc
      double complex stst2qq_M2Tree

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
	  
! Average factor
      if ((ftt1.eq.ftt2).and.(fgen1.eq.fgen2)) then
            avgfac = 1d0/(Nc**2)*1d0/2d0
      else 
            avgfac = 1d0/(Nc**2)
      end if
	 
      result(1) = intfac*flux*avgfac*stst2qq_M2Tree()

! ! !  ! ! ! ! ! For Resummation- singlet: 
      result(6) = intfac*flux*avgfac*stst2QQ_M2Tree3
! ! ! ! ! For Resummation- octet: 
      result(7) = intfac*flux*avgfac*stst2QQ_M2Tree6


	  end
	
*************************************

	  subroutine stst2qq_TreesplitDiffCS(result, costh)

	  implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac
      double precision result(13)
      double precision kappa  
      double complex stst2qq_M2Treesplit(13)!,try,stst2qq_M2Treesplittet
      external stst2qq_M2Treesplittet
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

! Average factor
      if ((ftt1.eq.ftt2).and.(fgen1.eq.fgen2)) then
      avgfac = 1d0/(Nc**2)*1d0/2d0
      else 
      avgfac = 1d0/(Nc**2)
      end if
      call stst2qq_M2Treesplittet(stst2qq_M2Treesplit)

      result(1) = intfac*flux*avgfac*stst2qq_M2Treesplit(1)
      result(2) = intfac*flux*avgfac*stst2qq_M2Treesplit(2)
      result(3) = intfac*flux*avgfac*stst2qq_M2Treesplit(3)
      result(4) = intfac*flux*avgfac*stst2qq_M2Treesplit(4)
      result(5) = intfac*flux*avgfac*stst2qq_M2Treesplit(5)
      result(6) = intfac*flux*avgfac*stst2qq_M2Treesplit(6)
      result(7) = intfac*flux*avgfac*stst2qq_M2Treesplit(7)
      result(8) = intfac*flux*avgfac*stst2qq_M2Treesplit(8)
      result(9) = intfac*flux*avgfac*stst2qq_M2Treesplit(9)
      result(10) = intfac*flux*avgfac*stst2qq_M2Treesplit(10)


      end

      subroutine stst2QQ_M2DipolePlusDist(intresdip, error)

      implicit none 
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
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
      double precision intresdip(ncomp), avgfac

      double precision RelAccuracy, AbsAccuracy,kappa
      double precision xi(ndim), fi(ncomp)

      external vegas 
      external stst2QQ_M2DipoleCS_Diff
 
 
 
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

      if (EpsPole.eq.0d0) then
 
 ! ===== VEGAS integration =====
      
      call vegas(ndim, ncomp, stst2QQ_M2DipoleCS_Diff,
     &      dble(RelAccuracy), dble(AbsAccuracy),
     &      flags, MinEval, MaxEval,
     &      Nstart, Nincrease,
     &      Neval, fail, intresdip, error, prob)
 

 ! ===== Average factor =====
! Flux factors for different units of cross-section
      if (fluxflag.eq.1) then
! in units of pb
cc	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,muu1**2*sqrtS**2,muu2**2*sqrtS**2)/(2d0*sqrtS))*sqrtS)
      ! in units of GeV^-2
      flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else
            
! in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
! 	 
      flux = 1d0/sqrtS**2
      endif
      
      ! Average factor
      if ((ftt1.eq.ftt2).and.(fgen1.eq.fgen2)) then
      avgfac = 1d0/(Nc**2)*1d0/2d0
      else 
      avgfac = 1d0/(Nc**2)
      end if
      intresdip=intresdip*flux*avgfac
      error=error*flux*avgfac
      endif
 
      end


      subroutine  stst2QQ_M2DipoleCS_Diff(ndim, xi, ncomp, fi)
      implicit none     
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"

      integer ndim, ncomp
      double precision xi(ndim), fi(ncomp) 
      double precision stst2QQ_M2DipolePlus_Quark,stst2QQ_M2DipolePlus_IEIS, intfac,sman

      sman = sqrtS**2

! Order of arguments ma2,mb2,mk2,mj2  
!******************* FE-IS/IE-FS: Integration for k1/p1 as emitter/spectator pair  ******************* 
      call stst2QQ_SetDipoleKinematicsFEIS(xi,x0FEIS,sman,m1**2,m2**2,m4**2,m3**2,intfac)

      fi(1) = intfac*stst2QQ_M2DipolePlus_Quark(1,1)

!******************* FE-IS/IE-FS: Integration for k1/p2 as emitter/spectator pair  ******************* 
      
      call stst2QQ_SetDipoleKinematicsFEIS(xi,x0FEIS,sman,m2**2,m1**2,m4**2,m3**2,intfac)

      fi(1) = fi(1) + intfac*stst2QQ_M2DipolePlus_Quark(2,0)

!******************* FE-IS/IE-FS: Integration for k2/p1 as emitter/spectator pair  ******************* 
      
      call stst2QQ_SetDipoleKinematicsFEIS(xi,x0FEIS,sman,m1**2,m2**2,m3**2,m4**2,intfac)

      fi(1) = fi(1) + intfac*stst2QQ_M2DipolePlus_Quark(2,0)

!******************* FE-IS/IE-FS: Integration for k2/p2 as emitter/spectator pair  ******************* 
      
      call stst2QQ_SetDipoleKinematicsFEIS(xi,x0FEIS,sman,m2**2,m1**2,m3**2,m4**2,intfac)

      fi(1) = fi(1) + intfac*stst2QQ_M2DipolePlus_Quark(1,1)

!******************* IE-FS:  p1 as emitter and p2 as spectator  ************************************************************
      call stst2QQ_SetDipoleKinematicsIEIS(xi,x0IEIS,sman,m1**2,m2**2,m3**2,m4**2,intfac) ! q^2=(p1-k1)^2=t

      fi(1) = fi(1) + intfac*stst2QQ_M2DipolePlus_IEIS()

!******************* IE-FS:  p2 as emitter and p1 as spectator  ************************************************************
      call stst2QQ_SetDipoleKinematicsIEIS(xi,x0IEIS,sman,m2**2,m1**2,m4**2,m3**2,intfac) ! q^2=(p2-k2)^2=t

      fi(1) = fi(1) + intfac*stst2QQ_M2DipolePlus_IEIS()

      end
c*************************************************************
	subroutine stst2QQ_Bremsstrahl23(intres, error)

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      parameter (ndim = 4)
      parameter (ncomp = 1)
      integer flags, mineval, maxeval, key, nregions, neval, fail
      integer maxiter, nstart, nincrease, nnew

      double precision flatness

      double precision intfac23, avgfac, x1min, x1max, x2min, x2max
      double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
      integer failcuh
      double precision error(ncomp), prob(ncomp)
      double precision intres(ncomp)
      double precision kappa

      double precision xgrid(0:50,ndim), xint, xacc(0:50,ndim), ymax(50,ndim), ymaxrat(50,ndim) , ans, err
      integer nhits(0:50,ndim), ifold(ndim), imode

      external integrand4
      external stst2QQ_integrand4
      external vegas
      external cuhre
      external suave
      external mint
      external stst2QQ_Integrand4Mint

	

      flags = 0
      mineval = 1000
      maxeval = 50000
      maxiter = 100
      nstart = maxeval/MAXITER
      nincrease = maxeval/(2*MAXITER)

      RELACCURACY = 5D-3
      ABSACCURACY = 1D-16

        
      call vegas(ndim, ncomp, stst2QQ_integrand4,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      flags, MINEVAL, MAXEVAL,
     &      nstart, nincrease,
     &      neval, fail, intres, error, prob)
        
c integration factor already includes the trivial integration over phi

      intfac23 = sqrtS**2/(32d0*(2*pi)**4)

c* flux 

      if (fluxflag.eq.1) then																!(Mo: new implemented by me (should be Patrick, because the comment was in german) , now processes 2->3 can be used without processes of 2->2, here are units of pb used (original:!Mo: neu von mir hinzugefuegt, jetzt kann auch der Prozess 2->3 ohne den Prozess 2->2 betrachtet werden, da hier )
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

      if ((ftt1.eq.ftt2).and.(fgen1.eq.fgen2)) then
      avgfac = 1d0/(Nc**2)*1d0/2d0
      else 
      avgfac = 1d0/(Nc**2)
      end if                                                              

      intres(1) = avgfac*flux*intfac23*intres(1)
      error(1) = avgfac*flux*intfac23*error(1)	

	end

!====================================================================
	subroutine stst2QQ_BremsstrahlWarnings(inptype, inpgen)
		
	implicit none
	
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
	
	integer inpgen, inptype
		
	!These warnings signal if the kinematics allow for a singular behavior of some t/u-channel propagator in the 2->3 part.In these cases, the integration routines often don't work correctly 
	!,which means that the estimated error is by far smaller than the actual uncertainty. Right now these cases are handeled by introducing a large width (see the Warning above).

! 	if(icha.eq.0) then
! 	
! !t,u
!            if((0.0.lt.(m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)).and.((m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)
!            print*,''
! 
!            elseif((0.0.lt.(m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)).and.((m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)
!            print*,''
! 
!            endif
! 
!           else
! !t
!            if((0.0.lt.(m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS)).and.((m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS)
!            print*,''
! 
!            elseif((0.0.lt.(m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS)).and.((m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS)
!            print*,''
! !u
!            elseif((0.0.lt.(m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)).and.((m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)
!            print*,''
! 
!            elseif((0.0.lt.(m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)).and.((m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS).lt.1.0)) then
! 
!            print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)
!            print*,''
! 
!            endif
! 
!           endif
		  
	end
			
	
	
!====================================================================
	subroutine stst2QQ_Integrand4(ndim, x, ncomp, f)

        implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
      integer ndim, ncomp
      double precision x(*), f(*)

      double precision eta,costh,x1max,x1min,x2max,x2min
      double precision cosxi, sinxi, sinth,stst2QQ_M2realemissionDipoles
      double complex stst2QQ_M2realemissionTGLTGL,stst2QQ_M2realemissionUGLUGL,stst2QQ_M2realemissionTGLUGL
      double complex stst2QQ_M2realemissionTNeuUNeu,stst2QQ_M2realemissionTNeuTNeu,stst2QQ_M2realemissionUNeuUNeu
      double complex stst2QQ_M2realemissionTGlUNeu,stst2QQ_M2realemissionUGlTNeu,stst2QQ_M2realemissionTGlTNeu,stst2QQ_M2realemissionUGlUNeu
      double complex stst2QQ_M2realemissionTGlTGl_prop, stst2QQ_M2realemissionUGlUGl_prop, stst2QQ_M2realemissionTGlUGl_prop, stst2QQ_M2realemissionTGlTNeu_prop
      double complex stst2QQ_M2realemissionUGlUNeu_prop, stst2QQ_M2realemissionTGlUNeu_prop, stst2QQ_M2realemissionUGlTNeu_prop


      call stst2QQ_PhSpace(x,muu3,muu4,muu5,muu1,muu2,cosxi)           
  
! ===== Single contributions to matrix element =====
      f(1) = 0d0
      if(dabs(cosxi).lt.0.999999999999999d0) then              !cut out numerical errors leading to a |cosxi|>1 oder cosxi=NAN (tests x1 and x2)

            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTGLTGL()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionUGLUGL()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTGLUGL()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTNeuTNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionUNeuUNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTNeuUNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTGlTNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionUGlUNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionTGlUNeu()
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionUGlTNeu()
    
            f(1) = f(1)+jacobi*stst2QQ_M2realemissionDipoles()
      
      endif

      end	
	
	
	
	
	
	
! ================= General Phase-space element======================
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine stst2QQ_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb,cosxi)

      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"

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
       	
       end
       
#include "stst2QQ_M2Tree.F"
#include "stst2QQ_M2Treesplittet.F"
#include "stst2QQ_M2Vertex.F"
#include "stst2QQ_M2Prop.F"
#include "stst2QQ_M2Counterterm.F"
#include "stst2QQ_M2CountertermsProp.F"
#include "stst2QQ_M2Boxitjs.F"
#include "stst2QQ_M2Boxjtis.F"
#include "stst2QQ_M2Boxtriplett.F"
#include "stst2QQ_M2Boxsextet.F"
#include "stst2QQ_M2Boxcancel3.F"
#include "stst2QQ_M2Boxcancel6.F"
#include "stst2QQ_M2realemissionTGlTGl.F"
#include "stst2QQ_M2realemissionUGlUGl.F"
#include "stst2QQ_M2realemissionTGlUGl.F"
#include "stst2QQ_M2realemissionTNeuUNeu.F"
#include "stst2QQ_M2realemissionTNeuTNeu.F"
#include "stst2QQ_M2realemissionUNeuUNeu.F"
#include "stst2QQ_M2realemissionTGlTNeu.F"
#include "stst2QQ_M2realemissionUGlUNeu.F"
#include "stst2QQ_M2realemissionTGlUNeu.F"
#include "stst2QQ_M2realemissionUGlTNeu.F"
#include "stst2QQ_M2Dipole.F"
#include "stst2QQ_M2realemissionDipoles.F"
#include "stst2QQ_Couplings.F"
#include "stst2QQ_Kinematics.F"
#include "stst2QQ_SetCounterterms.F"

