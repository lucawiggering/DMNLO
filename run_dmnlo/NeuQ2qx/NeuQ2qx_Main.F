
      subroutine neuq2qx(renscale,particle1,particle2,particle3,particle4,momentum,iflag,myresult,myerror,slhafilename)
      implicit none
#include "../util/DMNLO_LoopIntegrals.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
	
      integer ndim, ncomp, i,iflag(8),iflux, itree, imicr, ickm, MOswitch
      parameter (ncomp = 1)
      integer particle1,particle2,particle3,particle4
      double precision renscale,momentum
      double precision myresult(11),myerror(9), massth,result(9), error(9)
      double precision result23(ncomp),error23(ncomp)
      double precision intresdip(ncomp), errordip(ncomp),intresdip1(ncomp)
      character slhafilename*200

      QScale = renscale
      Qscalealphas = Qscale
 	
      Pcm = momentum
      n1 = particle1
      n2 = particle2
      n3 = particle3
      n4 = particle4
	
! ===== Flags assignment =====
      imass = iflag(1)
      iflux = iflag(2)
      itree = iflag(3)
      imicr = iflag(4)
      ickm = iflag(5) 
      ischeme = iflag(6)
      MOswitch = iflag(7)
      choosesol = iflag(8)

! ===== setting of masses with scale dependence =====	
      MTscale = 0
      MBscale = 0
	
! ===== Set model parameters =====
      call DMNLO_ModelPara(MOswitch,slhafilename)
      call DMNLO_ModelIni(slhafilename)   
! ===== Output model parameters =====
!      call DMNLO_ModelDigest

! ===== Renormalization Scales for Loop integrals =====
! EpsPoles = -2,-1,0,1,... steers LoopFunctions
! e.g. EpsPoles = -2/-1 lets all LoopFunctions return the IR diveregent coefficient corresponding to double/single pole
      EpsPole = 0d0
      UVdiv = 0d0
      IRdiv = 0d0	
        
* xiIR enables a check with mass regularization if set to 0d0 =====
* (default is xiIR = 1d0 !!!)
	xiIR = 1d0

! ===== Definition of Renormalization Scheme, Scale & Input =====
	call Init_RenScheme(ischeme)
      if(renfail.eq.0) stop
      
! ===== Kinematics and reduced masses =====
	call NeuQ2qx_SetKinematics()
* process energy
	sqrtS  = dsqrt(m1**2 + m2**2 + 2*Pcm**2 + 2*dsqrt(m1**2 + Pcm**2)*dsqrt(m2**2 + Pcm**2))
	sman = m1**2 + m2**2 + 2*Pcm**2 + 2*dsqrt(m1**2 + Pcm**2)*dsqrt(m2**2+ Pcm**2)

! Lower bound for x for FE-IS and IE-FS
	x0 =  (-m1**2 + m2**2 + m3**2 + (m1**2*m4)/sqrtS - (m4*(m2**2 - sman + m4*sqrtS))/sqrtS)/(-m1**2 + m2**2 + m4**2 + (m1**2*m4)/sqrtS + sman - 2*m4*sqrtS - (m4*(m2**2 - sman + m4*sqrtS))/sqrtS)

	if (x0.lt.0.9d0) then
		x0 =0.9d0
	endif

! The value of kappaDip can be chosen arbitrary. The final result must not depend on kappaDip.
	kappaDip=0d0

! ===== Check mass threshold =====
	massth = m3 + m4

	if (sqrtS.le.massth) then
      
      write(*,*) "NquQ2qx: overrun of mass threshold"

      do i=1,11
      myresult(i) = 0d0
      enddo   

      else 

      do i=1,9
      result(i) = 0d0
      enddo
      result23(1) = 0d0
	
* reduced masses
      m1red = m1/sqrtS
      m2red = m2/sqrtS
      m3red = m3/sqrtS
      m4red = m4/sqrtS
				
! ===== Phase-space slicing cut-off =====
!        cutS = 1e-6
!        cutC = 1e-5
!        if (n4.eq.21) then
!              cutoff = cutS*sqrtS
!              cutoffS = cutoff
!              cutoffC = cutC*sman
!        else if (n4.eq.22) then
!              cutoff = cutS*sqrtS
!              cutoffS = cutoff
!              cutoffC = 0d0
!        else
!              cutoff = cutS*sqrtS
!              cutoffS = 0d0
!              cutoffC = 0d0
!        endif
	
! ===== Set flux factor flag =====
      fluxflag = iflux
	
! ===== Generic couplings =====
!       (here the type and generation is set) 
      call NeuQ2qx_SetCouplings

! ===== Set Counterterms =====
      call NeuQ2qx_SetCounterterms

! ===== Integrate cross-section =====
      call NeuQ2qx_TreeCS(result, error, itree)

      call NeuQ2qx_M2DipoleCS(intresdip, errordip,itree)


      do i=1,9
      myerror(i) = error(i)
      enddo

      myerror(2) = errordip(1)

      result23(1) = 0d0

! 		Include bremsstrahlung
      if (itree.eq.0) call NeuQ2qx_Bremsstrahlung(result23, error23)

      do i=1,9
            myresult(i) = result(i)
      enddo
!    include x dep dipole part virtual contribution
      myresult(10) = result23(1)

      myresult(11) =  intresdip(1)

      myerror(1) = error23(1)

! add real part and plus distribution to the virtual part
      myresult(9) = myresult(9) + result23(1) + intresdip(1)
      endif

	end

! ========================================================================================================================================
	subroutine NeuQ2qx_TreeCS(result, error,mytree)

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"


	integer ncomp, ndim, i 
	
! ===== Set number of functions to integrate =====
      parameter (ncomp = 9)
      parameter (ndim = 1)

	double precision result(ncomp), error(ncomp), prob(ncomp)
	double precision costhmin, costhmax, relaccuracy, absaccuracy
	integer c, neval, fail, mytree
	integer flags, maxeval, mineval, maxiter, nstart, nincrease
	integer nnew, nregions, key
	integer nodes, npoints
	double precision flatness
	double precision myerror(ncomp)
	
!	double precision integrand, DGQUAD
	external NeuQ2qx_TreeDiffCS, NeuQ2qx_DiffCS, NeuQ2qx_IntegrandTree
	external Patterson,Gauss
	external suave, vegas, cuhre, divonne, DGQUAD, integrand
	
! ===== Set integration limits =====
	costhmin = -1d0
	costhmax = 1d0		
! ===== Settings =====

      flags = 0
      mineval = 1000
      maxeval = 100000
      maxiter = 30        
      nstart = maxeval/maxiter
      nincrease = maxeval/(2*maxiter)
      flatness = 50
      nnew = 1000
      key = 3
      nodes = 32
      npoints = 96
	
! ===== Set desired accuracy for integration =====
      absaccuracy = 1d-12
      relaccuracy = 1d-5
	
! ===== Compute full cross-section or only tree-level including channel contributions =====

	if (mytree.eq.0) then
	
!		Compute full cross-section
		call Patterson(ncomp, costhmin, costhmax, NeuQ2qx_DiffCS,
     &    		relaccuracy, absaccuracy, neval, fail, result, error)
     
!     		call vegas(ndim, ncomp, NeuQ2qx_IntegrandTree,
!     &                  dble(relaccuracy), dble(absaccuracy), flags, mineval, maxeval,
!     &                  nstart, nincrease, neval, fail, result, error, prob)
     
	else
	
!		Compute single tree-level contributions

		call Patterson(ncomp, costhmin, costhmax, NeuQ2qx_TreeDiffCS,
     &    		relaccuracy, absaccuracy, neval, fail, result, error)
     
     
!     		call vegas(ndim, ncomp, NeuQ2qx_IntegrandTree,
!     &                  dble(relaccuracy), dble(absaccuracy), flags, mineval, maxeval,
!     &                  nstart, nincrease, neval, fail, result, error, prob)

!		result(1) = DGQUAD(Integrand,costhmin,costhmax,npoints)


!   		call suave(ndim, ncomp, NeuQ2qx_IntegrandTree,
!     &                  dble(relaccuracy), dble(absaccuracy), flags, mineval, maxeval, nnew,
!     &                  flatness, nregions, neval, fail, result, error, prob)


!		call Gauss(ncomp, costhmin, costhmax, NeuQ2qx_TreeDiffCS, nodes, result)

!	  	call cuhre(ndim, ncomp, NeuQ2qx_IntegrandTree, dble(relaccuracy), dble(absaccuracy),
!     &      		 flags, mineval, maxeval, key, nregions, neval, fail, result, error, prob)
     
	endif

! ===== Check if accuracy has been reached =====
          
!	if( fail.ne.0 )
!     &    print *, "Failed to reach the desired accuracy."

	end

!*************************************************************************************

       subroutine NeuQ2qx_M2DipoleCS( intresdip, error,mytree)

       implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
     
! ***** Variables used for the integration method *****

       integer ndim, ncomp
       parameter (ndim = 2)
       parameter (ncomp = 1)
 

       integer flags, mineval, maxeval, key, nregions, neval, fail
       integer maxiter, nstart, nincrease,mytree

       double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
       integer failcuh
       double precision error(ncomp), prob(ncomp)
       double precision intresdip(ncomp), avgfac, avgcolor,kappa


       double precision RelAccuracy, AbsAccuracy
       double precision xi(ndim), fi(ncomp)

       external vegas 
       external NeuQ2qx_M2DipoleCS_DiffG,NeuQ2qx_M2DipoleCS_DiffEW



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

	   if (EpsPole.eq.0d0 .and. mytree.eq.0) then

! ===== VEGAS integration =====
               
	   if (xtype.eq.1) then

       call vegas(ndim, ncomp, NeuQ2qx_M2DipoleCS_DiffG,
     &      dble(RelAccuracy), dble(AbsAccuracy),
     &      flags, MinEval, MaxEval,
     &      Nstart, Nincrease,
     &      Neval, fail, intresdip, error, prob)

	   else 

		call vegas(ndim, ncomp, NeuQ2qx_M2DipoleCS_DiffEW,
     &      dble(RelAccuracy), dble(AbsAccuracy),
     &      flags, MinEval, MaxEval,
     &      Nstart, Nincrease,
     &      Neval, fail, intresdip, error, prob)

	   endif

! ===== Average factor =====
      avgfac = 1d0/2d0
      avgcolor = 1d0/NC	
	
! ===== Flux factor =====
 
 !	for sigma
	if (fluxflag.eq.1) then
!       	in units of GeV^-2
		flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
	else if (fluxflag.eq.2) then
!		in units of pb
		flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
	else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif

            intresdip=intresdip*flux*avgfac*avgcolor
            error = error*flux*avgfac*avgcolor 

      else

            intresdip = 0d0

      endif

       end

      subroutine NeuQ2qx_SetDipoleKinematics(xi,x0,sman,ma2,mmb2,mk2,mj2,intfac)
      implicit none 
! For the evaluation of the "plus"-distribution one has to distinguish between three different dependencies $f(x,Q^2)$. 
! For this purpose, we introduce two different variables for $Q^2$. These are q2 which is generated from the integration limits
! of $Q^2$ evaluated at x=1 $Q^2_\pm(1)$ and q2X which is generated from the $x$-dependent integration limits $Q^2_\pm(x)$ instead. 
! Variables with the suffix "X" refer to f(x,q2X), variables with the suffix "1" refer to f(1,q2) and those with the suffix "plus" correspond to f(x,q2).
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_DipoleKinematics.h"
      double precision xi(2),ma2,mmb2,mk2,mj2,intfac,kappa,alpha,beta,prefac,x0
      double precision q2max, q2min,q2maxAt1, q2minAt1,sman, zP, zM,rhoPlus, calAPlus,calBPlus

! set squared CM entergy for x=1
      sdip1 = sman

! generate value for x      
      xdip = (1d0-x0)*xi(1) + x0   

! integration limits for Q2
      prefac = 1d0/(2*(xdip*sman + (1d0-xdip)*(mmb2 - xdip*ma2)))

      alpha =  xdip**2*(ma2**2 + 2*ma2*(mmb2+mk2)-(mmb2-sman)**2)+ 2*mmb2*(ma2+mj2)
     &                 -xdip*(ma2**2 + ma2*(4*mmb2+mj2+mk2-sman) - (mmb2-sman)*(mmb2 - mj2 - mk2))

      beta =  xdip*kappa(ma2,mmb2,sman)*dsqrt((ma2-mmb2)**2*(1d0-xdip)**2 + (1d0-xdip)*(2*ma2*(mj2 + mk2*(2*xdip-1d0)-sman*xdip)
     &               -2*mmb2*(mj2+mk2-sman*xdip))+ kappa(xdip*sman,mk2,mj2)**2)

      q2max = prefac*(alpha + beta)
      q2min = prefac*(alpha - beta)

! q2minAt1 and q2maxAt1 correspond to the integration limits of Q2 evaluated at P2=mj2 which is equivalent to x=1
      q2minAt1 = mmb2 + mk2 - ((-ma2 + mmb2 + sman)*(mk2 - mj2 + sman))/(2*sman) - (kappa(sman,ma2,mmb2)*kappa(sman,mk2,mj2))/(2*sman)
      q2maxAt1 = mmb2 + mk2 - ((-ma2 + mmb2 + sman)*(mk2 - mj2 + sman))/(2*sman) + (kappa(sman,ma2,mmb2)*kappa(sman,mk2,mj2))/(2*sman)

! dQ2 is necessary for the substitution from q2 to q2X
      dQ2=(q2max-q2min)/(q2maxAt1-q2minAt1)

! q2 corresponds to Q2 for the integration limits evaluated at x=1. This variable has to be used in the part of the plus-distribution 
! which is proportional to the Dirac delta distribution
      q2 = (q2maxAt1-q2minAt1)*xi(2) + q2minAt1

! q2X is generated from the x-dependent integration limits
      q2X = dQ2*(q2-q2minAt1)+q2min
            
! jacobi factor 
      intfac = (q2maxAt1-q2minAt1)*(1d0-x0)/(8*Pi*kappa(sman,ma2,mmb2))

      if (q2X.gt.0d0) then
      if (xdip.le.(-q2X+ma2+mj2)/(2*dsqrt(ma2)*(dsqrt(ma2)-dsqrt(q2X)))) then
            write(*,*) "A higher value for x0 has to be chosen!"
      endif
      endif

! Define other variables that occur often and depend on x and Q2 
      q2bar1 = q2 - ma2 -mj2
      q2barX = q2X - ma2 - mj2

      sqrtLambda1 = kappa(q2,mj2,ma2)
      sqrtLambdaX = kappa(q2X,mj2,ma2)

      Rdip = dsqrt((q2barX + 2*ma2*xdip)**2 - 4*ma2*q2X*xdip**2)/sqrtLambdaX

      sdip = ma2 + mmb2 + 1d0/Rdip*(xdip*(sman-ma2-mmb2)+(q2barX+2*ma2*xdip)/(2*q2X)*(mmb2-mk2+q2X))
     &             - (q2X + ma2 -mj2)/(2*q2X)*(mmb2-mk2+q2X)

      etaA1 = ma2/(-q2bar1)
      etaAX = ma2/(-q2barX)

      etaJ1 = mj2/(-q2bar1)
      etaJX = mj2/(-q2barX)

      vDip1 = sqrtLambda1/(-q2bar1)
      vDipX = sqrtLambdaX/(-q2barX)

! Integration limits of z. These are useful for the definition of A(x) and B(x)
      zP = (1d0-xdip)/2d0*(-q2barX+ sqrtLambdaX*Rdip)/(xdip*mj2- q2barX*(1d0-xdip))
      zM = (1d0-xdip)/2d0*(-q2barX- sqrtLambdaX*Rdip)/(xdip*mj2- q2barX*(1d0-xdip))

! The following two variables are also global
      AdipX = (zP-zM)/(1d0-xdip+zM)
      BdipX = (zP-zM)/zM

      if (mj2.eq.0d0) then
      rhoPlus = dsqrt(1+ 4*etaA1*(xdip-1)*xdip)
      yAprimePlus = ((3 - 4*xdip)*etaA1 - 1)/rhoPlus**3
      calAPlus = 2*((1 - etaA1)*xdip-2)/(rhoPlus*(2*xdip - 3 - rhoPlus)) ! \mathcal{A}
      Aplus = 1d0/((1-xdip)*calAPlus)

      yBprimePlus = etaA1*(1 - 2*xdip)/rhoPlus**3
      calBPlus = 2*etaA1*xdip/(rhoPlus**2 + rhoPlus)  ! \mathcal{B}
      Bplus = 1d0/((1-xdip)*calBPlus)

      else
      Adip1 = 2*vdip1/(2*etaJ1+1-vdip1)
      Bdip1 = 2*vdip1/(1-vdip1)
      endif

      end
      
!
!=========================================================================================================================================

      subroutine  NeuQ2qx_M2DipoleCS_DiffEW(ndim, xi, ncomp, fi)
      implicit none 
 
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      integer ndim, ncomp
      double precision xi(ndim), fi(ncomp) 
      double precision NeuQ2qx_M2DipolePlusEW, intfac

      call NeuQ2qx_SetDipoleKinematics(xi,x0,sman,m2**2,m1**2,m4**2,m3**2,intfac)
	  
      fi(1) = intfac*NeuQ2qx_M2DipolePlusEW()
!
      end

      subroutine  NeuQ2qx_M2DipoleCS_DiffG(ndim, xi, ncomp, fi)
      implicit none 

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      integer ndim, ncomp
      double precision xi(ndim), fi(ncomp) 
      double precision NeuQ2qx_M2DipolePlusG_Quark,NeuQ2qx_M2DipolePlusG_Gluon,intfac

! quark as emitter
      call NeuQ2qx_SetDipoleKinematics(xi,x0,sman,m2**2,m1**2,m4**2,m3**2,intfac)
      
      fi(1) = intfac*NeuQ2qx_M2DipolePlusG_Quark()

! gluon as emitter
      call NeuQ2qx_SetDipoleKinematics(xi,x0,sman,m2**2,m1**2,m3**2,m4**2,intfac)     

      fi(1) = fi(1)+intfac*NeuQ2qx_M2DipolePlusG_Gluon()

      end
! ============================================================================================================================================



	subroutine NeuQ2qx_DiffCS(result, costh)
	
	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
	double precision costh, intfac, intfacMO, avgfac, avgcolor,costheta, result(9), kappa
	double precision NeuQ2qx_M2Tree,  NeuQ2qx_M2Box, NeuQ2qx_M2Vertex, NeuQ2qx_M2Propagator
	double precision NeuQ2qx_M2VertexCounter, NeuQ2qx_M2PropagatorCounter,NeuQ2qx_M2DipoleEW,NeuQ2qx_M2DipoleG, NeuQ2qx_M2TreePerc
	double precision NeuQ2qx_SoftPhoton, NeuQ2qx_SoftGluon, NeuQ2qx_HardCollGluon

	double precision NeuQ2qx_M2TreeDipole

      result = 0d0  
! ===== Set Mandelstam variables t and u =====

      tman = (m1**2 + m3**2 - (sqrtS**2 + m3**2 - m4**2)*(sqrtS**2 + m1**2 - m2**2)/(2*sqrtS**2)+costh*kappa(sqrtS**2, m3**2, m4**2)*kappa(sqrtS**2, m1**2, m2**2)/(2*sqrtS**2))
      tred = (m1red**2 + m3red**2 - (1 + m3red**2 - m4red**2)*(1 + m1red**2- m2red**2)/2+costh*kappa(1d0, m3red**2, m4red**2)*kappa(1d0, m1red**2, m2red**2)/2)
      uman = (m1**2 + m2**2 + m3**2 + m4**2 - sman - tman)
      ured = (m1red**2 + m2red**2 + m3red**2 + m4red**2 - 1d0 - tred)
      
! ===== Integration factor =====
!       (the azimuthal integration constant 2pi is inlcuded)

	intfac = 2*pi*(kappa(sqrtS**2,m3**2,m4**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
! ===== Flux factor =====
!
!	for sigma
      if (fluxflag.eq.1) then
!       	in units of GeV^-2
            flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
      else if (fluxflag.eq.2) then
!		in units of pb
            flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
      else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
      flux = 1d0/sqrtS**2
      endif

! ===== Average factor =====

	avgfac = 1d0/2d0
	
! ===== Average factor =====

	avgcolor = 1d0/NC	
	
! ===== Single contributions to the cross-section =====
!       (no flux and averaging factor)
!	Tree-level
	result(1) = 		intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()

! 	Full NLO contribution
	result(2) = 0d0
!        if(cutoff.lt.zerocheck) then
      if (xtype .eq. 1) then ! gluon in the final state
            result(2) =            intfac*flux*avgfac*avgcolor*NeuQ2qx_M2DipoleG()
      else !  EW particles in the final state
            result(2) =            intfac*flux*avgfac*avgcolor*NeuQ2qx_M2DipoleEW()
      endif 
!        else 
!        result(2) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()*NeuQ2qx_SoftGluon()
!        endif
     
! 	Soft photon divergence substraction (only for photon final state)
!      if (n4.eq.22) then
!            n4 = 21
!            call NeuQ2qx_SetKinematics()
!            result(2) = result(2) + intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()*NeuQ2qx_SoftPhoton()
!            n4 = 22
!            call NeuQ2qx_SetKinematics()
!      endif

!	Vertex contributions
	result(3) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Vertex()
	
!	Propagator contributions
	result(4) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Propagator()
!    model of hard coll virtual part 2*B(0,0,0)
! 	result(4) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()*((GS**2*CV)/(16*pi**2))*2*(-1d0)

!	Vertex counterterms
	result(5) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2VertexCounter()

!	Propagator counterterms
	result(6) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2PropagatorCounter()

!	Boxen result
	result(7) = 	        intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Box()
	
!	Hard collinear part of gluon final state	
	result(8) = 0d0
!	result(8) =             intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()*NeuQ2qx_HardCollGluon()
	
!	sum
	result(9) = sum(result(1:8))
      end

! =====================================================================================================================================
	subroutine NeuQ2qx_TreeDiffCS(result, costh)

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision costh, intfac, intfacMO, avgfac, avgcolor
	double precision result(9)
	double precision kappa, NeuQ2qx_M2Tree,  NeuQ2qx_M2TreeH, NeuQ2qx_M2TreeV, NeuQ2qx_M2TreePerc, NeuQ2qx_M2Vertex, NeuQ2qx_M2VertexH, NeuQ2qx_M2VertexV, NeuQ2qx_M2PropagatorCounter
	double precision NeuQ2qx_M2Propagator, NeuQ2qx_M2VertexCounter
	
! ===== Mandelstam variables t,u =====

	tman = (m1**2 + m3**2 - (sqrtS**2 + m3**2 - m4**2)*(sqrtS**2 + m1**2 - m2**2)/(2*sqrtS**2)+costh*kappa(sqrtS**2, m3**2, m4**2)*kappa(sqrtS**2, m1**2, m2**2)/(2*sqrtS**2))
	tred = ((m1red**2 + m3red**2 - (1 + m3red**2 - m4red**2)*(1 + m1red**2 - m2red**2)/2d0)+(costh*kappa(1d0, m3red**2, m4red**2)*kappa(1d0, m1red**2, m2red**2)/2d0))
	uman = (m1**2 + m2**2 + m3**2 + m4**2 - sman - tman)
	ured = (m1red**2 + m2red**2 + m3red**2 + m4red**2 - 1d0 - tred)
	
! ====== Prefactors ======
!        (the azimuthal integration constant 2pi is included in intfac)

	intfac = 2*pi*(kappa(sman,m3**2,m4**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	intfacMO = 2*pi*(kappa(sman**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	
! ====== Flux factor =====

!	for sigma
      if (fluxflag.eq.1) then
!       	in units of GeV^-2
            flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
      else if (fluxflag.eq.2) then
!		in units of pb
            flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
      else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif

! ===== Average factor =====

      avgfac = 1d0/2d0
                  
! ===== Average factor =====

      avgcolor = 1d0/NC
			
! ===== Results including flux and averaging factor =====

      result(1) =   intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()
      result(2) =   intfac*flux*avgfac*avgcolor*treeperc(1)
      result(3) =   intfac*flux*avgfac*avgcolor*treeperc(2)
      result(4) =   intfac*flux*avgfac*avgcolor*treeperc(3)
      result(5) =   intfac*flux*avgfac*avgcolor*treeperc(4)
      result(6) =   intfac*flux*avgfac*avgcolor*treeperc(5)
      result(7) =   intfac*flux*avgfac*avgcolor*treeperc(6)
      result(8) =   intfac*flux*avgfac*avgcolor*treeperc(7)
      result(9) = 0
      
	end
! ====================================================================================================================================
	subroutine NeuQ2qx_IntegrandTree(ndim, x, ncomp, f)
	
	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	integer ndim, ncomp
	double precision x(ndim), f(ncomp)
	double precision x1max,x1min
	double precision costh, intfac, intfacMO, avgfac, avgcolor
	double precision kappa, NeuQ2qx_M2Tree, NeuQ2qx_M2TreePerc
	
! ===== Integration boundaries and Jacobi determinant =====

	x1min = -1d0
	x1max = 1d0
	costh = (x1max-x1min)*x(1)+x1min
	jacobi = x1max-x1min
	

! ===== Mandelstam variables t,u =====

	tman = (m1**2 + m3**2 - (sqrtS**2 + m3**2 - m4**2)*(sqrtS**2 + m1**2 - m2**2)/(2*sqrtS**2)+costh*kappa(sqrtS**2, m3**2, m4**2)*kappa(sqrtS**2, m1**2, m2**2)/(2*sqrtS**2))
	tred = ((m1red**2 + m3red**2 - (1 + m3red**2 - m4red**2)*(1 + m1red**2 - m2red**2)/2d0)+(costh*kappa(1d0, m3red**2, m4red**2)*kappa(1d0, m1red**2, m2red**2)/2d0))
	uman = (m1**2 + m2**2 + m3**2 + m4**2 - sman - tman)
	ured = (m1red**2 + m2red**2 + m3red**2 + m4red**2 - 1.0 - tred)
	
! ====== Prefactors ======
!        (the azimuthal integration constant 2pi is included in intfac)

	intfac = 2*pi*(kappa(sman,m3**2,m4**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	intfacMO = 2*pi*(kappa(sman**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	
! ====== Flux factor =====

!	for sigma
      if (fluxflag.eq.1) then
!       	in units of GeV^-2
            flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
      else if (fluxflag.eq.2) then
!		in units of pb
            flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
      else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif

! ===== Average factor =====

	avgfac = 1d0/2d0
		
	
! ===== Average factor =====

	avgcolor = 1d0/NC	
		
! ===== Results including flux and averaging factor =====

	f(1) = 0d0
	f(1) = jacobi*intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()
	f(2) = 0d0
	f(3) = 0d0
	f(4) = 0d0
	f(5) = 0d0
	f(6) = 0d0
	f(7) = 0d0
	f(8) = 0d0
    
	end	
			
! =========================================================================================================================================
      subroutine NeuQ2qx_Bremsstrahlung(intres, error)

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      integer ndim, ncomp
      parameter (ndim = 4)
      parameter (ncomp = 1)
      integer flags, mineval, maxeval, key, nregions, neval, fail
      integer maxiter, nstart, nincrease

      double precision intfac23,x1min, x1max, x2min, x2max, kappa
      double precision intcuh(ncomp), ercuh(ncomp), probcuh(ncomp)
      integer failcuh
      double precision error(ncomp), prob(ncomp)
      double precision intres(ncomp)
      double precision RelAccuracy, AbsAccuracy
      double precision x(ndim)




      external NeuQ2qx_Integrand,NeuQ2qx_Integrand_Test
      external vegas
      external cuhre
	
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

! ===== CUHRE integration =====
!	  call cuhre(ndim, ncomp, NeuQ2qx_Integrand_Test,
!     &      dble(RelAccuracy), dble(AbsAccuracy),
!     &      flags, MinEval, MaxEval,
!     &      Nstart, Nincrease,
!     &      Neval, fail, intres, error, prob)

!	intcuh(1) = intres(1)
!	ercuh(1) = error(1)
!	failcuh = fail
!	probcuh(1) = prob(1)

! ===== VEGAS integration =====

 
        
	call vegas(ndim, ncomp, NeuQ2qx_Integrand_Test,
     &      dble(RelAccuracy), dble(AbsAccuracy),
     &      flags, MinEval, MaxEval,
     &      Nstart, Nincrease,
     &      Neval, fail, intres, error, prob)
	

! ===== Flux factor =====

!	for sigma
      if (fluxflag.eq.1) then
!       	in units of GeV^-2
            flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
      else if (fluxflag.eq.2) then
!		in units of pb
            flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
      else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif


! ===== Integration factor =====
!       (already includes the trivial integration over phi)

	intfac23 = sqrtS**2/(32d0*(2*pi)**4)


! ===== Final cross-section =====

	intres(1) = flux*intfac23*intres(1)
	error(1)  = flux*intfac23*error(1)


	end

! ====================================================================================================================================
	subroutine NeuQ2qx_Integrand_Test(ndim, x, ncomp, f)

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
	integer ndim, ncomp
	double precision x(*), f(*)

	double precision eta,costh,x1max,x1min,x2max,x2min
	double precision ifac, neutfac,avgcolor,avgfac
	double precision cosxi, sinxi, sinth
	double precision NeuQ2qx_M2Radiation, NeuQ2qx_M2RadiationGqq
	double precision NeuQ2qx_M2RadiationGSTcon, NeuQ2qx_M2RadiationGTTcon, NeuQ2qx_M2RadiationDipoleEW
	double precision NeuQ2qx_M2RadiationGSSsoft, NeuQ2qx_M2RadiationGSTsoft, NeuQ2qx_M2RadiationGTTsoft
	double precision NeuQ2qx_M2RadiationGSSsoftcoll, NeuQ2qx_M2RadiationGSTsoftcoll, NeuQ2qx_M2RadiationGTTsoftcoll,NeuQ2qx_M2RadiationDipoleG
	double precision NeuQ2qx_M2RadiationGSScoll, NeuQ2qx_M2RadiationGSTcoll, NeuQ2qx_M2RadiationGTTcoll,NeuQ2qx_M2RadiationDipoleGqq
	
	external NeuQ2qx_PhSpace

! ===== Integration boundaries and Jacobi determinant =====

! dipole configuration -> gluon integrated out
!      cutoff = 0d0
!      call NeuQ2qx_PhSpace(x,m3red,0d0,m4red,m1red,m2red)

! soft-gluon configuration -> no gluon integrated out
! cutoff defined in NeuQ2qx_Main.F line 81
       
        call NeuQ2qx_PhSpace(x,m3red,m4red,0d0,m1red,m2red)


! ===== Scalar products - dipole configuration =====
!         x1 = xx1
!         x2 = xx3
       
! 	pap1 = kak1
! 	pap2 = kak3             
! 	pbp1 = kbk1
! 	pbp2 = kbk3

! ===== Scalar products - soft-gluon configuration =====
      x1 = xx1
      x2 = xx2


      papb = kakb
      pap1 = kak1
      pap2 = kak2  
      pap3 = kak3

      pbp1 = kbk1
      pbp2 = kbk2
      pbp3 = kbk3

      p1p2 = k1k2
      p1p3 = k1k3
      p2p3 = k2k3


! ===== Averaging factor =====

	avgfac = 1d0/2d0

! ===== Average factor =====

	avgcolor = 1d0/NC

! ===== Single contributions to matrix element =====

        f(1) =  0d0
!
! if gluon
!
        if (xtype.eq.1) then  ! gluon
! convergent contribution          
           f(1) = jacobi*avgfac*avgcolor*(NeuQ2qx_M2RadiationGSTcon() + NeuQ2qx_M2RadiationGTTcon())
! pure soft contribution     
!           if ((xx2.gt.(2d0*cutoffS/dsqrt(sman))).and.(xx3.gt.(2d0*cutoff/dsqrt(sman)))) then
            f(1) = f(1) + jacobi*avgfac*avgcolor*(NeuQ2qx_M2RadiationGSSsoft() + NeuQ2qx_M2RadiationGSTsoft() + NeuQ2qx_M2RadiationGTTsoft())
!           endif
! softcollinear contribution
!           if ((xx1.lt.(1d0 + m3red**2 - cutoffC/sman)).and.((xx2.gt.(2d0*cutoffS/dsqrt(sman))).and.(xx3.gt.(2d0*cutoff/dsqrt(sman))))) then
           f(1) = f(1) + jacobi*avgfac*avgcolor*(NeuQ2qx_M2RadiationGSSsoftcoll() + NeuQ2qx_M2RadiationGSTsoftcoll() + NeuQ2qx_M2RadiationGTTsoftcoll())
!           endif
! collinear contribution
!           if (xx1.lt.(1d0 + m3red**2 - cutoffC/sman)) then
            f(1) = f(1) + jacobi*avgfac*avgcolor*(NeuQ2qx_M2RadiationGSScoll() + NeuQ2qx_M2RadiationGSTcoll() + NeuQ2qx_M2RadiationGTTcoll())
		f(1) = f(1) + jacobi*avgfac*avgcolor*(NeuQ2qx_M2RadiationGqq() + NeuQ2qx_M2RadiationDipoleGqq())
!           endif
! dipole gluon
		f(1) = f(1) + jacobi*avgfac*avgcolor*NeuQ2qx_M2RadiationDipoleG()

!           
! all other particles
! < only soft divergences >
!       
      else 
      
            f(1) = jacobi*avgfac*avgcolor*(NeuQ2qx_M2Radiation() + NeuQ2qx_M2RadiationDipoleEW())

      endif
           

	end

! ================= General Phase-space element
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine NeuQ2qx_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb)

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      double precision x(*)
      double precision mmr1,mmr2,mmr3,mmra,mmrb
      double precision eta, costh
      double precision cosxi,sinxi,sinth
      double precision xx3max,xx3min,xx1max,xx1min
      double precision taubar,sigmabar,muplus,muminus
      

! ===== Integration boundaries and Jacobi determinant (from notes) =====

      xx3min = 2*mmr3
                  
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
      
      if (cosxi.lt.(-1d0)) then
            cosxi = -1d0
      endif
      
      if (cosxi.gt.(1d0)) then
            cosxi = 1d0
      endif

      sinxi = dsqrt(1-cosxi**2)
      sinth = dsqrt(1-costh**2)
      
      
! ===== Scalar products =====

      kakb = (sman/2d0)*(1d0-mmra**2-mmrb**2)
      
      k1k2 = (sman/2d0)*(1d0 - xx3 - mmr1**2 - mmr2**2 + mmr3**2)
      k1k3 = (sman/2d0)*(1d0 - xx2 - mmr1**2 + mmr2**2 - mmr3**2)
      k2k3 = (sman/2d0)*(1d0 - xx1 + mmr1**2 - mmr2**2 - mmr3**2)

      kak1 = (sman/4d0)*(xx1*(1+mmra**2-mmrb**2)-dsqrt((1d0-mmra**2-mmrb**2)**2-4*mmra**2*mmrb**2)*dsqrt(xx1**2-4*mmr1**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
      kak3 = (sman/4d0)*(xx3*(1+mmra**2-mmrb**2)-dsqrt((1d0-mmra**2-mmrb**2)**2-4*mmra**2*mmrb**2)*dsqrt(xx3**2-4*mmr3**2)*costh) 
      kak2 = (sman/2d0)*(1d0+mmra**2-mmrb**2)-kak1-kak3
      
      kbk1 = (sman/4d0)*(xx1*(1-mmra**2+mmrb**2)+dsqrt((1d0-mmra**2-mmrb**2)**2-4*mmra**2*mmrb**2)*dsqrt(xx1**2-4*mmr1**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
      kbk3 = (sman/4d0)*(xx3*(1-mmra**2+mmrb**2)+dsqrt((1d0-mmra**2-mmrb**2)**2-4*mmra**2*mmrb**2)*dsqrt(xx3**2-4*mmr3**2)*costh)
      kbk2 = (sman/2d0)*(1d0-mmra**2+mmrb**2)-kbk1-kbk3	


      end


! ===== Integrand Function for Tree Level =====

	double precision function integrand(costh)

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision costh, intfac, intfacMO, avgfac, avgcolor
	double precision kappa, NeuQ2qx_M2Tree,  NeuQ2qx_M2TreePerc

! ===== Mandelstam variables t,u =====

	tman = (m1**2 + m3**2 - (sqrtS**2 + m3**2 - m4**2)*(sqrtS**2 + m1**2 - m2**2)/(2*sqrtS**2)+costh*kappa(sqrtS**2, m3**2, m4**2)*kappa(sqrtS**2, m1**2, m2**2)/(2*sqrtS**2))
	tred = ((m1red**2 + m3red**2 - (1 + m3red**2 - m4red**2)*(1 + m1red**2 - m2red**2)/2d0)+(costh*kappa(1d0, m3red**2, m4red**2)*kappa(1d0, m1red**2, m2red**2)/2d0))
	uman = (m1**2 + m2**2 + m3**2 + m4**2 - sman - tman)
	ured = (m1red**2 + m2red**2 + m3red**2 + m4red**2 - 1.0 - tred)

! ====== Prefactors ======
!        (the azimuthal integration constant 2pi is included in intfac)

	intfac = 2*pi*(kappa(sman,m3**2,m4**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	intfacMO = 2*pi*(kappa(sman**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)
	
! ====== Flux factor =====

!	for sigma
      if (fluxflag.eq.1) then
!       	in units of GeV^-2
            flux = 1d0/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
      else if (fluxflag.eq.2) then
!		in units of pb
            flux = hbar_c2/(4d0*(kappa(sman,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
!	for sigma*v
      else
! 		in units of cm^3/sec
!	 	flux = 2.9979d-26*hbar_c2/sqrtS**2
! 		in units of GeV^-2
            flux = 1d0/sqrtS**2
      endif

! ===== Average factor =====

	avgfac = 1d0/2d0
		
	
! ===== Average factor =====

	avgcolor = 1d0/NC
	
		
! ===== Results including flux and averaging factor =====
	integrand = intfac*flux*avgfac*avgcolor*NeuQ2qx_M2Tree()
	

	end
		
! =======================================================================================================================================
! Matrix elements to include (in seperate files) follow here:

#include "NeuQ2qx_Couplings.F"
#include "NeuQ2qx_Kinematics.F"
	

#include "NeuQ2qx_M2Tree.F"
#include "NeuQ2qx_M2TreeGA.F"
#include "NeuQ2qx_M2TreeGAred.F"
#include "NeuQ2qx_M2TreeH.F"
#include "NeuQ2qx_M2TreeHgen.F"
#include "NeuQ2qx_M2TreeV.F"
#include "NeuQ2qx_M2TreeVgen.F"
#include "NeuQ2qx_M2TreePerc.F"


#include "NeuQ2qx_M2AmpH.F"
#include "NeuQ2qx_M2AmpBoxH.F"

#include "NeuQ2qx_M2Vertex.F"
#include "NeuQ2qx_M2VertexH.F"
#include "NeuQ2qx_M2VertexV.F"
#include "NeuQ2qx_M2VertexG.F"
#include "NeuQ2qx_M2Propagator.F"
#include "NeuQ2qx_M2PropagatorH.F"
#include "NeuQ2qx_M2PropagatorV.F"
#include "NeuQ2qx_M2Box.F"
#include "NeuQ2qx_M2BoxH.F"
#include "NeuQ2qx_M2BoxV.F"
#include "NeuQ2qx_M2BoxG.F"
#include "NeuQ2qx_M2Radiation.F"
#include "NeuQ2qx_M2RadiationH.F"
#include "NeuQ2qx_M2RadiationV.F"
#include "NeuQ2qx_M2RadiationG.F"
#include "NeuQ2qx_M2RadiationDipole.F"
!#include "NeuQ2qx_SoftGluon.F"

#include "NeuQ2qx_SetCounterterms.F"
!#include "NeuQ2qx_SetCountertermsV.F"
#include "NeuQ2qx_M2Counterterms.F"
#include "NeuQ2qx_M2PropagatorCounter.F"
#include "NeuQ2qx_M2PropagatorCounterH.F"
#include "NeuQ2qx_M2PropagatorCounterV.F"
#include "NeuQ2qx_M2VertexCounter.F"
#include "NeuQ2qx_M2VertexCounterH.F"
#include "NeuQ2qx_M2VertexCounterV.F"
#include "NeuQ2qx_M2VertexCounterG.F"

#include "NeuQ2qx_M2DipoleInt.F"
#include "NeuQ2qx_M2Dipole.F"

