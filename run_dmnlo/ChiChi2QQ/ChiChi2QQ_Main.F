      subroutine chichi2qq(renscale,Pcm,iflag,cut,myresult,slhafilename)
      implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_GenCouplings.h" 
#include "ChiChi2QQ_GenCounterterms.h"

      integer ndim, ncomp, iflag(14), inptype, inpgen, iflux, i, j, k
      parameter (ncomp = 1)

      double precision Pcm,myresult(9), sinth, costh, avgfac, massth
      double precision renscale
      double precision result(9), error(9)
      double precision result23(ncomp),error23(ncomp),cut
      character slhafilename*200
      integer MOswitch
      

      QScale = renscale
      Qscalealphas = Qscale

** flags assignment
c icharg = 0- only neutralinos in initial state, 1- chargino & neutralino mixed, 2 - charginos only (chargino+ is p_1, chargino- is p_2)
c icha = number of charginos in initial state (0...2)
	icha = iflag(1)

! indi = index of neutralino or chargino (1...4 for neutralino, 1...2 for chargino)
	indi = iflag(2)
	
! indj = index of the neutralino (1...4 for neutralino, 1...2 for chargino, typically indj=1)
	indj = iflag(3)
      
! inptype = 1 - neutrino, 2- lepton, 3- up quark, 4- down quark (in chargino neutralino co-ann. type is of the up-type lepton/quark)
	inptype = iflag(4)

	if(icha==1.AND.(inptype==4.OR.inptype==2)) then
		write(*,*)'WARNING: In Chargino Neutralino Co-Annihilation the final state has to be up-type!'
	endif
	
! inpgen = generation of the quarks
	inpgen = iflag(5)

! imass = 0- take masses from MicrOmegas, 1-diagonalize here	
	imass = iflag(6)

! iflux = 0 - calculate v.sigma, 1- standard cross-section sigma	
	iflux = iflag(7)

! itree = 0- full one loop, 1- only tree level (quicker)	
	itree = iflag(8)

! imicr = kinematics in MicrOmegas cross-section - 0-our best estimate, 1- MO original 	
	imicr = iflag(9)
	
! Patrick: To choose between the two options for the Resummation, see also ChiChi2QQ_Counterterms.f. 
! (0 = Old DMNLO variant, 1 = Newer Spira variant)
	Res_Style = iflag(10)
	
! Patrick: Switch the NNLO Resummation on and off (1 = on, 0 = off). The use NNLO_Res = 1 sets Res_Style to 1 automatically.
	NNLO_Res = iflag(11)

! ischeme - choose renormalization scheme
      ischeme = iflag(12)

! legacy
      MOswitch = iflag(13)

! choosesol
      choosesol = iflag(14)
	

! ===== Set model parameters =====
! reading the MSSM parameters from micrOmegas and setting the masses & mixings & all shortforms e.g. TB,SA
	call DMNLO_ModelPara(MOswitch,slhafilename)
	call DMNLO_ModelIni(slhafilename)
!	call DMNLO_ModelDigest

! ===== set sfermion-widths to zero/10 ===========
      do i = 1,2 
      do j = 3,4
      do k = 1,3

      if(WSf(i,j,k).lt.10.0) then
      !see Warning below for purpose of this
      !print*,'WARNING: MicrOMEGAs sQuark-width small (=',WSf(i,j,k),'GeV )=> set width WSf(',i,',',j,',',k,') to 10 GeV'
      WSf(i,j,k)=10.0
      endif

      enddo
      enddo
      enddo

! Define the CMS \sqrt{s} and the scale Q2
	if (icha.eq.0) then
		sqrtS  = dsqrt(Pcm**2 + MNeu2(indi)) + dsqrt(Pcm**2 + MNeu2(indj))
	else if (icha.eq.1) then
		sqrtS  = dsqrt(Pcm**2 + MCha2(indi)) + dsqrt(Pcm**2 + MNeu2(indj))
	else if (icha.eq.2) then
		sqrtS  = dsqrt(Pcm**2 + MCha2(indi)) + dsqrt(Pcm**2 + MCha2(indj))
	endif
! QScale old value used in ntnt2bb, now QScale = dsqrt(MSf(1,3,3)*MSf(2,3,3)) set in Init_RenScheme 	
!        QScale = sqrtS        
	
! ===== Renormalization Scales for Loop integrals =====
! Divergent UV & IR poles for Loop integrals
      UVdiv = 0d5
      IRdiv = 0d3

! ===== Dipole (cutoff=0) or soft Gluon (cutoff >10^(-10)) =====
      cutoff = cut*sqrtS
	
! Renormalization scale \mu & switch between dim.reg and mass regularization for IR divergence - xiIR
! xiIR enables a check with mass regularization if set to 0d0 -> DEFAULT is xiIR = 1d0 !!!! 
      xiIR = 1d0
      muSc = QScale

! EpsPoles = -2,-1,0,1,... steers LoopFunctions
! e.g. EpsPoles = -2/-1 lets all LoopFunctions return the IR diveregent coefficient corresponding to double/single pole
! checking UV divergence requires setting EpsPole=1 and UVdiv=1 (loop amplitude should be zero or extremely small)
      EpsPole = 0d0

! ===== Definition of Renormalization Scheme, Scale & Input =====
!      Set which solution should be taken for the Mixmatrices USf - only implemented in Renscheme 1 (0:DMNLO chooses solution, 1: take solution1, 2: take solution 2)
      call Init_RenScheme(ischeme)
      if(renfail.eq.0) stop
      
! set the fermion type of the final state and flux factor (in mykinematics.h)
! in case of chargino-neutralino annihilation tt = type of the final state *anti*particle 
      tt = inptype
      gen = inpgen
      fluxflag = iflux
      
! Generic couplings (here the type and generation is set !!! )
	call ChiChi2QQ_SetCouplings
	
! FA kinematics & reduced masses
	call ChiChi2QQ_SetKinematics
	
! self-energies & counterterms
	call ChiChi2QQ_SetCounterTerms()
	
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
	call ChiChi2QQ_IntCS(result,error,itree)

! Compute differential cross-section
cc call ChiChi2QQ_DiffCS(result, dcos(pi/4d0))
	result23(1) = 0d0

! Compute bremsstrahlung contribution       
      if(itree.eq.0) then

      call ChiChi2QQ_BremsstrahlWarnings(inptype, inpgen)
      call ChiChi2QQ_Bremsstrahl23(result23, error23) 

      endif
	
! Copy to output vectors
         
      if (itree.eq.0) then

      do i=1,6
      myresult(i) = result(i)
      enddo

      myresult(7) = result23(1)
      myresult(8) = error23(1)
      myresult(9) = result(6) + result23(1)

      else

      do i=1,9
      myresult(i) = result(i)
      enddo

      endif

      endif
	
	end

**************************************************************

	subroutine ChiChi2QQ_IntCS(result,error,mytree)

	implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"

	double precision result(9), error(9)

	double precision costhmin, costhmax, relaccuracy, absaccuracy
	double precision kappa
	integer c, neval, fail, ncomp, mytree
	external ChiChi2QQ_TreeDiffCS, ChiChi2QQ_DiffCS, Patterson

        ncomp = 7

! Integration limits for cos(theta)
	costhmin = -1d0
	costhmax = 1d0

! Set requested accuracy
	absaccuracy = 1d-12
	relaccuracy = 1d-5

! Full calculation
	if (mytree.eq.0) then

	call Patterson(ncomp, costhmin, costhmax, ChiChi2QQ_DiffCS, 
     &    relaccuracy, absaccuracy, neval, fail, result, error)

! Only tree-level calculation
	else

	call Patterson(ncomp, costhmin, costhmax, ChiChi2QQ_TreeDiffCS,
     &    relaccuracy, absaccuracy, neval, fail, result, error)

	endif

c	if( fail .ne. 0 )
c     &    print *, "Failed to reach the desired accuracy."
	
	end

*************************************

	subroutine ChiChi2QQ_DiffCS(result, costh)

	implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	double precision costh, intfac, intfacMO, avgfac
	double precision result(9)

	double precision kappa
	double precision ChiChi2QQ_M2Tree, ChiChi2QQ_M2TreeM0, ChiChi2QQ_M2TreeHO, ChiChi2QQ_M2BoxGluino, ChiChi2QQ_M2BoxGluon
	double precision ChiChi2QQ_M2Vertex, ChiChi2QQ_M2Propagator, ChiChi2QQ_SoftGluon
	double precision ChiChi2QQ_M2VertexSCounter, ChiChi2QQ_M2VertexTUCounter, ChiChi2QQ_M2PropagatorCounter, ChiChi2QQ_M2Dipole

	double precision test

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
	intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

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
      avgfac = 1d0/4d0

! ***************************************
! contributions to the matrix element 
! ***************************************
! (1-tree, 2-full one loop virtual, 3-vertex, 4-propagator, 5-box, 6-dipole, 7-micrOmegas tree, 8-higher order effective)

      result(1) = intfac*flux*avgfac*ChiChi2QQ_M2Tree()

      result(2) = intfac*flux*avgfac*ChiChi2QQ_M2Vertex()
      result(2) = result(2) + intfac*flux*avgfac*ChiChi2QQ_M2VertexSCounter()
      result(2) = result(2) + intfac*flux*avgfac*ChiChi2QQ_M2VertexTUCounter()

      result(3) = intfac*flux*avgfac*ChiChi2QQ_M2Propagator()
      result(3) = result(3) + intfac*flux*avgfac*ChiChi2QQ_M2PropagatorCounter()

      result(4)  = intfac*flux*avgfac*ChiChi2QQ_M2BoxGluon()
      result(4) = result(4) + intfac*flux*avgfac*ChiChi2QQ_M2BoxGluino()	  
	  
      if(cutoff.lt.zerocheck) then

      result(5) = intfac*flux*avgfac*ChiChi2QQ_M2Dipole()

      else  

      result(5) = intfac*flux*avgfac*ChiChi2QQ_SoftGluon()*ChiChi2QQ_M2Tree()

      endif         

      result(6) = result(1) + result(2) + result(3) + result(4) + result(5) + intfac*flux*avgfac*ChiChi2QQ_M2TreeHO()
	  
	end

*************************************

	subroutine ChiChi2QQ_TreeDiffCS(result, costh)

	implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	double precision costh, intfac, intfacMO, avgfac
	double precision result(9)

	double precision kappa, ChiChi2QQ_M2Tree, ChiChi2QQ_M2TreeM0, ChiChi2QQ_M2TreeHO, ChiChi2QQ_M2TreePerc

! Mandelstam variables t=(p2-k1)**2 , u=(p1-k1)**2
! Scattering angle theta between p1 and k1
c CHECK CONVENTION T/U ????

	tman = sqrtS**2*(muu2**2 + muu3**2 - 0.5d0*(1d0-muu1**2+muu2**2)*(1d0-muu4**2+muu3**2) - 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
	tred =          (muu2**2 + muu3**2 - 0.5d0*(1d0-muu1**2+muu2**2)*(1d0-muu4**2+muu3**2) - 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

	uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)
	ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

! Integration factor for 2->2 phase-space integration in CMS 
! the azimuthal integration constant 2pi is inlcuded in intfac
	intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)
	intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

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
      avgfac = 1d0/4d0
	
! no flux and averaging factor!

      result(1) = intfac*flux*avgfac*ChiChi2QQ_M2Tree()
      result(2) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(11)  !V*V
      result(3) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(12)  !V*S
      result(4) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(22)  !S*S
      result(5) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(13)  !V*t/u
      result(6) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(23)  !S*t/u
      result(7) = intfac*flux*avgfac*ChiChi2QQ_M2TreePerc(33)  !t/u*t/u

	end


c*************************************************************
	subroutine ChiChi2QQ_Bremsstrahl23(intres, error)

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

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
      external ChiChi2QQ_integrand4
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

	  call vegas(ndim, ncomp, ChiChi2QQ_integrand4,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      flags, MINEVAL, MAXEVAL,
     &      nstart, nincrease,
     &      neval, fail, intres, error, prob)


        
c integration factor already includes the trivial integration over phi

	intfac23 = sqrtS**2/(32d0*(2*pi)**4)

c* flux 

      if (fluxflag.eq.1) then																!Mo: neu von mir hinzugefuegt, jetzt kann auch der Prozess 2->3 ohne den Prozess 2->2 betrachtet werden, da hier 
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

	avgfac = 1d0/4d0                                                                

	intres(1) = avgfac*flux*intfac23*intres(1)
	error(1) = avgfac*flux*intfac23*error(1)

	end

!====================================================================
	subroutine ChiChi2QQ_BremsstrahlWarnings(inptype, inpgen)
		
	implicit none
	
#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
	
	integer inpgen, inptype
		
	!These warnings signal if the kinematics allow for a singular behavior of some t/u-channel propagator in the 2->3 part.In these cases, the integration routines often don't work correctly 
	!,which means that the estimated error is by far smaller than the actual uncertainty. Right now these cases are handeled by introducing a large width (see the Warning above).

	if(icha.eq.0) then
	
!t,u
           if((0.0.lt.(m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)).and.((m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m4**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)
           print*,''

           elseif((0.0.lt.(m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)).and.((m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m3**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)
           print*,''

           endif

          else
!t
           if((0.0.lt.(m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS)).and.((m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m4**2-MSf2(1,inptype+1,inpgen))/(m2*sqrtS)
           print*,''

           elseif((0.0.lt.(m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS)).and.((m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m3**2-MSf2(1,inptype+1,inpgen))/(m1*sqrtS)
           print*,''
!u
           elseif((0.0.lt.(m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)).and.((m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m2**2+m3**2-MSf2(1,inptype,inpgen))/(m2*sqrtS)
           print*,''

           elseif((0.0.lt.(m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)).and.((m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS).lt.1.0)) then

           print*,'WARNING: Kinematics may lead to an unstable 2->3-integration!',(m1**2+m4**2-MSf2(1,inptype,inpgen))/(m1*sqrtS)
           print*,''

           endif

          endif
		  
	end
		
		

!====================================================================
	subroutine ChiChi2QQ_Integrand4(ndim, x, ncomp, f)

        implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
	integer ndim, ncomp
	double precision x(*), f(*)

	double precision eta,costh,x1max,x1min,x2max,x2min
	double precision cosxi, sinxi, sinth
	double precision ChiChi2QQ_RadiationS,ChiChi2QQ_RadiationT,ChiChi2QQ_RadiationU,ChiChi2QQ_RadiationTS,ChiChi2QQ_RadiationUS,ChiChi2QQ_RadiationTU
	
	external ChiChi2QQ_PhSpace

! ===== Integration boundaries and Jacobi determinant =====

! dipole configuration -> gluon integrated out
      if(cutoff.lt.zerocheck) then


       call ChiChi2QQ_PhSpace(x,muu3,muu5,muu4,muu1,muu2,cosxi)            !exchange masses as needed for correct 2->3 phase-space integration; integration over Quark1 and Quark2 phase space (for quark-quark-Gluon-finalstate)

! ===== Scalar products - dipole configuration =====

      x1 = xx1
      x2 = xx3    

      else 

!       print*,'softgluon'

! soft-gluon configuration -> no gluon integrated out
! cutoff defined in ChiChi2QQ_Main.F line 87
       call ChiChi2QQ_PhSpace(x,muu3,muu4,muu5,muu1,muu2,cosxi)            !exchange masses as needed for correct 2->3 phase-space integration; integration over Quark1 and Gluon phase space (for quark-quark-Gluon-finalstate)

! ===== Scalar products - soft-gluon configuration =====

      x1 = xx1
      x2 = xx2

      endif
	
! ===== Single contributions to matrix element =====

      if(dabs(cosxi).lt.0.999999999999999d0) then              !cut out numerical errors leading to a |cosxi|>1 oder cosxi=NAN (tests x1 and x2)

	f(1) = 0d0					    !Die Wirkungsquerschnitte in Radiation-Datei enthalten integriert Anteile
                                                            !=> Ausgabe des differentiellen WQ's also nicht mehr möglich !!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      f(1) = jacobi*ChiChi2QQ_RadiationS()   
      f(1) = f(1) + jacobi*ChiChi2QQ_RadiationT()
      f(1) = f(1) + jacobi*ChiChi2QQ_RadiationU()
      f(1) = f(1) + jacobi*ChiChi2QQ_RadiationTS()
      f(1) = f(1) + jacobi*ChiChi2QQ_RadiationUS()               
      f(1) = f(1) + jacobi*ChiChi2QQ_RadiationTU()    

      else

      f(1) = 0d0

      endif


	end

! ================= General Phase-space element======================
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine ChiChi2QQ_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb,cosxi)

      implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"

      double precision x(*)
      double precision mmr1,mmr2,mmr3,mmra,mmrb
      double precision eta, costh,p1k3
      double precision cosxi,sinxi,sinth
      double precision xx3max,xx3min,xx1max,xx1min
      double precision taubar,sigmabar,muplus,muminus


! ***** Integration boundaries and Jacobi determinant (from HitchHiker) *****

      if (cutoff.lt.zerocheck) then
            xx3min = 2*mmr3              !k3->k2
      else 
            xx3min = 2d0*cutoff/sqrtS    !k3->k3,also einfach eine reduzierte Masse des Gluons als cutoff eingeführt
      endif
                  
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

       if(cutoff.lt.zerocheck) then

!Dipol subtraction

       p1k1=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
       p1k2=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)   !Quark-Gluon-Skalarprodukt (aber p1k2 ist im Hauptcode 
                                                                                                                                           !als QuarkQuark-Skalarprodukt verwendet=> p1k2<->p1k3
               
       else

!Phase space slicing
      
       p1k1=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
       p1k3=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)  !Quark-Gluon-Skalarprodukt (aber p1k2 ist im Hauptcode als QuarkQuark-Skalarprodukt verwendet
                                                                                                                                          !=>umschreiben notwendig!!
       p1k2=(sqrtS**2/2d0)*(1d0-mmrb**2+mmra**2)-p1k1-p1k3                                                                                !Impulserhaltung: p1k2=p1p1+p1p2-p1k1-p1k3

       endif

!******************************************************************************


	end

#include "ChiChi2QQ_Couplings.F"
#include "ChiChi2QQ_Kinematics.F"
#include "ChiChi2QQ_M2Tree.F"
#include "ChiChi2QQ_M2TreeMO.F"
#include "ChiChi2QQ_M2TreeHO.F"
#include "ChiChi2QQ_M2TreePerc.F"
#include "ChiChi2QQ_M2Dipole.F"
#include "ChiChi2QQ_SoftGluon.F"
#include "ChiChi2QQ_M2Vertex.F"
#include "ChiChi2QQ_M2Propagator.F"
#include "ChiChi2QQ_M2GluonBox.F"
#include "ChiChi2QQ_M2GluinoBox.F"
#include "ChiChi2QQ_M2Counterterms.F"
#include "ChiChi2QQ_M2Radiation.F"           
#include "ChiChi2QQ_SetCounterterms.F"
