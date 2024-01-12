* QQ2xx.F
* sQuark (co-) annihilations into elw. final states 
      subroutine qq2xx(renscale,particle1,particle2,particle3,particle4,particle5,momentum,iflag,myresult,myerror,slhafilename)
      implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"   


      integer iflag(8),inptype, inpgen,itree, iflux, imicr, i,MOswitch
      integer particle1,particle2,particle3,particle4,particle5,io_error
      double precision renscale,myresult(10), myerror(8), massth,cut, momentum                                                                               
      double precision result(8), error(8)
      double precision result23(1), error23(1) , result_Coulomb
      double complex QQ2xx_M2CoulombCorrection
      character slhafilename*200



** flags assignment

! imass = 0- take masses from MicrOmegas, 1-diagonalize here	
      imass = iflag(1)

! iflux = 0 - calculate v.sigma, 1- standard cross-section sigma	
      iflux = iflag(2)

! itree = 0- full one loop, 1- only tree level (quicker)	
      itree = iflag(3)

! imicr = kinematics in MicrOmegas cross-section - 0-our best estimate, 1- MO original 	
      imicr = iflag(4)

! ischeme = renormalizatino scheme
      ischeme = iflag(5)

! isomm = 1 - add Sommerfled enhancement, 0 don't
      isomm = iflag(6) ! in Model.h

! legacy 
      MOswitch = iflag(7) ! in Model.h

! choosesol 
      choosesol = iflag(8)


!set fluxflag
      fluxflag = iflux

! scales
      QScale = renscale
      Qscalealphas = Qscale
	
!set particle-types and Momentum
 
      pcm = momentum
      n1 = particle1
      n2 = particle2
      n3 = particle3
      n4 = particle4
      n5 = particle5
         
!************************************************************************************************************************         
         
!choose channel and corrections        

      channel = 99    !full tree level

      corrProp = 99   !full Propagator correction
      corrVert = 99   !full Vertex correction
      corrBox = 99    !full Box correction
         
! for debugging (to switch on velocity enhanced terms only in vertex and box)

!          corrProp = 999  !Sommerfeld correction
!          corrVert = 999  !Sommerfeld correction
!          corrBox = 999   !Sommerfeld correction

!*************************************************************************************************************************

! reading the MSSM parameters from micrOmegas and setting the masses & mixings & all shortforms e.g. TB,SA

!old                                 (must be taken + switch off Init_RenScheme(ischeme) for tree-level-comparison with CalcHEP only!)
!	call QQ2xx_ModelPara
!	call QQ2xx_ModelIni(imass)

!new
      call DMNLO_ModelPara(MOswitch,slhafilename)
      call DMNLO_ModelIni(slhafilename)
!      call DMNLO_ModelDigest

!*************************************************************************************************************************

! Divergent UV & IR poles for Loop integrals

      UVdiv = 0d0
      IRdiv = 0d0    
	
! Renormalization scale \mu & switch between dim.reg and mass regularization for IR divergence - xiIR
! xiIR enables a check with mass regularization if set to 0d0 -> DEFAULT is xiIR = 1d0 !!!!
 
      xiIR = 1d0
!	muSc = QScale

! EpsPoles = -2,-1,0,1,... steers LoopFunctions
! e.g. EpsPoles = -2/-1 lets all LoopFunctions return the IR diveregent coefficient corresponding to double/single pole
! checking UV divergence requires setting EpsPole=1 and UVdiv=1 (loop amplitude should be zero or extremely small)
      EpsPole = 0d0

! Definition of Renormalization Scheme, Scale & Input
      call Init_RenScheme(ischeme)        
      if(renfail.eq.0) stop                            

! FA kinematics & masses of external particles
       
	call QQ2xx_SetKinematics(itree)
     
! Define the CMS \sqrt{s} and the scale Q2

	sqrtS  = dsqrt(Pcm**2 + m1**2) + dsqrt(Pcm**2 + m2**2)
        
!Phase-space slicing cut-off 
      cut = 0.000005
	cutoff = cut*sqrtS

!set reduced masses

      muu1=m1/sqrtS
      muu2=m2/sqrtS
      muu3=m3/sqrtS
      muu4=m4/sqrtS
      muu5=m5/sqrtS

! Generic couplings (here the type and generation is set !!! )

      call QQ2xx_SetCouplings                                  

! set counterterms 

	call QQ2xx_SetCounterTerms

*******************************************

! Mass threshold
      massth = m3+m4
      if (sqrtS.le.massth) then

            do i=1,8
            result(i) = 0d0
            enddo

            result_Coulomb = 0d0
            result23(1) = 0d0

            print*,'final state not accessable'

      else 

! 2->2 cross-section
! flag abort22: 3-terminate without calculation (result=0), 
!               2-on-shell subtraction case, 
!               1-full 2->2 calculation

          if(abort22.lt.3) then

! Compute integrated cross-section
            
            if(abort22.eq.0) then

	      call QQ2xx_IntCS(result,error,itree)   !full-1-loop (itree=0) or tree level (itree=1) calculation

              result_Coulomb = 0d0

              if((itree.eq.0).and.(isomm.eq.1)) result_Coulomb = QQ2xx_M2CoulombCorrection(0,0) * result(1)  !CoulombCorrections(*,*) first argument: 0 = standard method, 1=approximation via improved Sommerfeld-Sakharov factor                                                     !                        second argument: 0= NLO-subtraction from NLO-Coulomb, 1= NLO-subtraction from full NLO calculation
                                                                                          !                        if corrProp = corrVert = corrBox = 999-> (*,1) for NLO cross-check 
              if((itree.eq.1).and.(isomm.eq.1)) result_Coulomb = QQ2xx_M2CoulombCorrection(0,2) * result(1)
            else

	      call QQ2xx_IntCS(result,error,1)       !tree-level-calculations with (abort22=2) or without (abort22=1) on-shell-subtraction  

              result_Coulomb = 0d0

            endif         
       
! Compute differential cross-section
!         call QQ2xx_TreeDiffCS(result, 0.5d0)  !Check of differential cross sections at tree level
!	 call QQ2xx_DiffCS(result, dcos(pi/3d0))
c      print*,'test result=',result(1)

          else

	    do i=1,8
	      result(i) = 0d0
	    enddo

            result_Coulomb = 0d0

          endif
!2->3
          if(abort23.eq.0) then

	    result23(1) = 0d0
            error23(1)  = 0d0

            if (itree.eq.0) call QQ2xx_Bremsstrahl23(result23, error23)  

          else

            result23(1) = 0d0      
            error23(1)  = 0d0

          endif


        endif

! Copy to output vectors
	do i=1,6
            myresult(i) = result(i)
            myerror(i) = error(i)
	enddo

!Coulomb
      myresult(7)  = result_Coulomb
!        myresult(8)  = result_Coulomb/5d0     !very rough estimate of the error due to scale unceratinties/unknown higher order corrections (taken at the Coulomb-peak) (see arXiv:0812.0919v2 p.15, 0501289v2 p.17 )

!2->3
      myresult(9)  = result23(1)
      myresult(10) = error23(1)

!NLO
      myresult(8)= myresult(2) + result23(1)
      myerror(8) = dsqrt(myerror(2)**2 + error23(1)**2)                       

!Full result
      myresult(2) = myresult(2) + result23(1) + result_Coulomb   
      myerror(2) = dsqrt(myerror(2)**2 + error23(1)**2 + myresult(7)**2/25d0)

	end

!*******************************************************************************************

	subroutine QQ2xx_IntCS(result,error,mytree)

	implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

      integer ndim, ncomp

      parameter (ndim =4)    !8 !Taken form NeuQ2qx
      parameter (ncomp = 8)   

      double precision result(ncomp), error(ncomp), prob(ncomp)

      double precision costhmin, costhmax, relaccuracy, absaccuracy
      integer c, neval, fail, mytree
      integer flags, maxeval, mineval, maxiter, nstart, nincrease
      integer nnew, nregions, key
      integer nodes, npoints
      double precision flatness
      external QQ2xx_TreeDiffCS, QQ2xx_DiffCS, QQ2xx_IntegrandTree
      external Patterson, Gauss
      external vegas, suave, cuhre, DGQUAD, integrand

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
	absaccuracy = 1d-12!1d-8
	relaccuracy = 1d-5!1d-2

! Full calculation

	if (mytree.eq.0) then

	call Patterson(ncomp, costhmin, costhmax, QQ2xx_DiffCS,
     &    relaccuracy, absaccuracy, neval, fail, result, error)

! Only tree-level calculation
	else

	call Patterson(ncomp, costhmin, costhmax, QQ2xx_TreeDiffCS,   
     &    relaccuracy, absaccuracy, neval, fail, result, error)

!    	call vegas(ndim, ncomp, QQ2xx_IntegrandTree,
!     &                  dble(relaccuracy), dble(absaccuracy), flags, mineval, maxeval,
!     &                  nstart, nincrease, neval, fail, result, error, prob)
           
!	result(1) = DGQUAD(Integrand,costhmin,costhmax,npoints)


!   	call suave(ndim, ncomp, QQ2xx_IntegrandTree,
!     &                  dble(relaccuracy), dble(absaccuracy), flags, mineval, maxeval, nnew,
!     &                  flatness, nregions, neval, fail, result, error, prob)


!	call Gauss(ncomp, costhmin, costhmax, QQ2xx_TreeDiffCS, nodes, result)

!	call cuhre(ndim, ncomp, QQ2xx_IntegrandTree, dble(relaccuracy), dble(absaccuracy),
!     &  flags, mineval, maxeval, key, nregions, neval, fail, result, error, prob)   

      endif

!	if( fail .ne. 0 )
!     &    print *, "Failed to reach the desired accuracy."

      end

*************************************

      subroutine QQ2xx_DiffCS(result, costh)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac, S, color
      double precision result(8), QQ2xx_M2Tree, QQ2xx_M2Prop, QQ2xx_M2Vert, QQ2xx_M2Box, QQ2xx_SoftGluon, QQ2xx_M2Dipole
      double precision kappa

      external QQ2xx_M2AmpInit2to2

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
!	intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

! Flux factors for different units of cross-section
	if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
            flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)
      else

! in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
	 flux = 1d0/sqrtS**2
	endif

! Average factor

        if((xtype.eq.ytype).and.(xtype.ne.0).and.(xtype.ne.4).and.(xtype.ne.8)) then

!Necessary, because phase-space-integration includes the points (p1,q2) and (p2,q1) (1,2 are the particles , p and q
!two momenta) which double counts the same process if particle 1 and 2 are identical.
          S=2d0                                                                               

                                                                                             
        elseif((f3type.eq.f4type).and.(f3type.ne.0).and.(f3gen.eq.f4gen).and.(total_chargeout.ne.0)) then
              
!covers the case of quark-quark or antiquark-antiquark final state
         S=2d0                                                                               

        else

         S=1d0

        endif

! Color Average factor
        
        if(((sftype1.eq.3).or.(sftype1.eq.4)).and.((sftype2.eq.3).or.(sftype2.eq.4))) then
         
         Color = NC
           
           if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then
            
!divide by 9 (1/3 more because elw. final-state-particles are not colored,
!but fermions are multiplied by NC in the corresponding amplitude)
            Color = Color*NC                                                                 
                                                                                             
           endif 

!covers case of non colored particles in the final state AND in the initial state
        else if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then   
         
!divide by 3 (1/3 more because elw. final-state-particles are not colored,
!but are multiplied by NC in the corresponding amplitude)
         Color = NC                                                                          
                                                                                             
        else
         
         Color=1d0
        
        endif

	avgfac = 1d0/(S*Color)                                                               !only scalars incoming, S=Symmetryfactor(=2 if final-state-particles are the same, 
	
! ***************************************
! contributions to the matrix element 
! ***************************************
! (1-tree, 2-full one loop virtual, 3-vertex, 4-propagator, 5-box, 6-dipole/phase-space-slicing)

        call QQ2xx_M2AmpInit2to2(0,channel,corrProp,corrVert,corrBox,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)

!******************************************************

	result(1) = intfac*flux*avgfac*QQ2xx_M2Tree

!******************************************************

        result(2) = 0d0 

        result(2) = result(1)                                                                                                                                                                        

        result(2) = result(2) + intfac*flux*avgfac*QQ2xx_M2Vert       

        result(2) = result(2) + intfac*flux*avgfac*QQ2xx_M2Prop

        result(2) = result(2) + intfac*flux*avgfac*QQ2xx_M2Box

c  	 result(2) = result(2) + intfac*flux*avgfac*QQ2xx_M2Dipole()

        result(2) = result(2) + QQ2xx_SoftGluon()*result(1) 

!******************************************************	

	result(3) = intfac*flux*avgfac*QQ2xx_M2Vert                      

!	print*,'Mvert=',QQ2xx_M2Vert

!******************************************************

	result(4) = intfac*flux*avgfac*QQ2xx_M2Prop                 

!	print*,'Mprop=',QQ2xx_M2Prop

!******************************************************
c
	result(5) = intfac*flux*avgfac*QQ2xx_M2Box

!	print*,'Mbox=',QQ2xx_M2Box

!******************************************************

!Dipoles
!        result(6) = intfac*flux*avgfac*QQ2xx_M2Dipole()

!Soft Gluon (phase-space-slicing)
        result(6)= QQ2xx_SoftGluon()*result(1)

	end

!******************************************************

	subroutine QQ2xx_TreeDiffCS(result, costh)

	implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"

      double precision costh, intfac, intfacMO, avgfac, Color, S
      double precision result(8)

      double precision kappa, QQ2xx_M2Tree, QQ2xx_M2Prop, QQ2xx_M2Vert, QQ2xx_M2Box, QQ2xx_M2TreePerc

      external QQ2xx_M2AmpInit2to2

! Mandelstam variables t=(p1-k1)**2 , u=(p1-k2)**2
! Scattering angle theta between p1 and k1

	tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) +   
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
	tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

	uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)                       
	ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

! Integration factor for 2->2 phase-space integration in CMS 
! the azimuthal integration constant 2pi is inlcuded in intfac
	intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)
!	intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

! Flux factor
	if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
      flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)

      else

! in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
	 flux = 1d0/sqrtS**2
	endif

        if((xtype.eq.ytype).and.(xtype.ne.0).and.(xtype.ne.4).and.(xtype.ne.8)) then

         S=2d0                                                                               !Necessary, because phase-space-integration includes the points (p1,q2) and (p2,q1) (1,2 are the particles , p and q
                                                                                             !two momenta) which doublecounts the same process if particle 1 and 2 are identical.
        elseif((f3type.eq.f4type).and.(f3type.ne.0).and.(f3gen.eq.f4gen).and.(total_chargeout.ne.0)) then

         S=2d0                                                                               !covers the case of quark-quark or antiquark-antiquark final state

        else

         S=1d0

        endif

        if(((sftype1.eq.3).or.(sftype1.eq.4)).and.((sftype2.eq.3).or.(sftype2.eq.4))) then
         Color = NC
           if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then
            Color = Color*NC                                                                 !non-colored final state: divide by 9 (1/3 more because elw. final-state-particles are not colored, but amplitudes multiplied by NC=3)
           endif 
        elseif(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then
         Color = NC                                                                          !non-colored final AND initial state: divide by 3 (1/3 more because elw. final-state-particles are not colored, but amplitudes multiplied by NC=3)
        else
         Color=1d0
        endif

	avgfac = 1d0/(S*Color)                                                               !only scalars incoming=> no spin-averaging
                                                        
!only usable for tree-level-checks

!           call QQ2xx_M2AmpInit2to2(1,1,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !Q*Q
! 	  result(1) = intfac*flux*avgfac*QQ2xx_M2Tree
! !           print*,'QQ',result(1)
! 
!           call QQ2xx_M2AmpInit2to2(1,2,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !Q*S
! 	  result(2) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
! !           print*,'QS',result(2)
! 
!           call QQ2xx_M2AmpInit2to2(1,3,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !Q*V
! 	  result(3) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'QV',result(3)
! 
!           call QQ2xx_M2AmpInit2to2(1,4,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !Q*T
!   	  result(4) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
!           call QQ2xx_M2AmpInit2to2(1,5,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !Q*U
!   	  result(4) = result(4) + 2*intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'Q*(T+U)',result(4)
! 
!           call QQ2xx_M2AmpInit2to2(1,7,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !S*S
! 	  result(5) = intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'SS',result(5)
! 
!          call QQ2xx_M2AmpInit2to2(1,99,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !full tree
! 	  result(6) = intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'full tree',result(6)

 
!            call QQ2xx_M2AmpInit2to2(1,8,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !S*V
!            result(1) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
!  !          print*,'SV',result(1)
! 
!           call QQ2xx_M2AmpInit2to2(1,9,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)    !S*T
! 	  result(2) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
!           call QQ2xx_M2AmpInit2to2(1,10,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !S*U
! 	  result(2) = result(2) + 2*intfac*flux*avgfac*QQ2xx_M2Tree
! !           print*,'S*(T+U)',result(2)
! 
!           call QQ2xx_M2AmpInit2to2(1,13,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !V*V
! 	  result(3) = intfac*flux*avgfac*QQ2xx_M2Tree
! !           print*,'VV',result(3)
! 
!           call QQ2xx_M2AmpInit2to2(1,14,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !V*T
! 	  result(4) = 2*intfac*flux*avgfac*QQ2xx_M2Tree
!           call QQ2xx_M2AmpInit2to2(1,15,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !V*U
! 	  result(4) = result(4) + 2*intfac*flux*avgfac*QQ2xx_M2Tree
! !           print*,'V*(T+U)',result(4)
! 
!           call QQ2xx_M2AmpInit2to2(1,19,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !T*T
!   	  result(5) = intfac*flux*avgfac*QQ2xx_M2Tree
!           call QQ2xx_M2AmpInit2to2(1,20,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !T*U
!   	  result(5) = result(5) + 2*intfac*flux*avgfac*QQ2xx_M2Tree
!           call QQ2xx_M2AmpInit2to2(1,25,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !U*U
!   	  result(5) = result(5) + intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'(T+U)*(T+U)',result(5)
! 
!           call QQ2xx_M2AmpInit2to2(1,99,99,99,99,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)   !full tree
! 	  result(6) = intfac*flux*avgfac*QQ2xx_M2Tree
! !          print*,'full tree',result(6)

!default choice

      call QQ2xx_M2AmpInit2to2(1,channel,corrProp,corrVert,corrBox,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)

      result(1) = intfac*flux*avgfac*QQ2xx_M2Tree

      result(2) = intfac*flux*avgfac*QQ2xx_M2Tree      !tree in result(2) is needed for DMNLOscan.cpp

      result(3) = 0d0
      result(4) = 0d0
      result(5) = 0d0
      result(6) = 0d0
      result(7) = 0d0
      result(8) = 0d0

	end


c*************************************************************

	subroutine QQ2xx_IntegrandTree(ndim, x, ncomp, f)

	implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"


        integer ndim, ncomp
	double precision x(ndim), f(ncomp)      !ndim=1, ncomp=8
	double precision x1max,x1min
	double precision costh, intfac, intfacMO, avgfac, Color, S

	double precision kappa, QQ2xx_M2Tree, QQ2xx_M2Prop, QQ2xx_M2Vert, QQ2xx_M2Box

        external QQ2xx_M2AmpInit2to2

! Mandelstam variables t=(p1-k1)**2 , u=(p1-k2)**2
! Scattering angle theta between p1 and k1

! ***** Integration boundaries and Jacobi determinant *****

	x1min = - 1d0
	x1max = 1d0
	costh = (x1max-x1min)*x(1)+x1min
	jacobi = x1max-x1min

!***********************************************************

	tman = sqrtS**2*(muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) +  
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)
	tred =          (muu1**2 + muu3**2 - 0.5d0*(1d0+muu1**2-muu2**2)*(1d0-muu4**2+muu3**2) + 
     &	           0.5d0*kappa(1d0,muu1**2,muu2**2)*kappa(1d0,muu3**2,muu4**2)*costh)

	uman = sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred)                       
	ured =           muu1**2 + muu2**2 + muu3**2 + muu4**2 - 1d0 - tred

! Integration factor for 2->2 phase-space integration in CMS 
! the azimuthal integration constant 2pi is inlcuded in intfac

	intfac = 2*pi*(kappa(1d0,muu3**2,muu4**2)/2d0)/(4*(2*pi)**2)
!	intfacMO = 2*pi*(kappa(sqrtS**2,Mfm(4)**2,Mfm(4)**2)/(2d0*sqrtS))/(4*(2*pi)**2*sqrtS)

! Flux factor
	if (fluxflag.eq.1) then
! in units of pb
!	 flux = hbar_c2/(4d0*(kappa(sqrtS**2,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
! in units of GeV^-2
	 flux = 1d0/(4d0*(kappa(1d0,muu1**2,muu2**2)/2d0)*sqrtS**2)

      else

! in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
! in units of GeV^-2
	 flux = 1d0/sqrtS**2
	endif

        if((xtype.eq.ytype).and.(xtype.ne.0).and.(xtype.ne.4).and.(xtype.ne.8)) then

         S=2d0                                                                               !Necessary, because phase-space-integration includes the points (p1,q2) and (p2,q1) (1,2 are the particles , p and q
                                                                                             !two momenta) which doublecounts the same process if particle 1 and 2 are identical.
        elseif((f3type.eq.f4type).and.(f3type.ne.0).and.(f3gen.eq.f4gen).and.(total_chargeout.ne.0)) then

         S=2d0                                                                               !covers the case of quark-quark or antiquark-antiquark final state
                                                                       
        else

         S=1d0

        endif

        if(((sftype1.eq.3).or.(sftype1.eq.4)).and.((sftype2.eq.3).or.(sftype2.eq.4))) then
         Color = NC
           if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then
            Color = Color*NC                                                                 !non-colored final state: divide by 9 (1/3 more because elw. final-state-particles are not colored, but amplitudes multiplied by NC=3)
           endif 
        elseif(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then
         Color = NC                                                                          !non-colored final AND initial state: divide by 3 (1/3 more because elw. final-state-particles are not colored, but amplitudes multiplied by NC=3)
        else
         Color=1d0
        endif

	avgfac = 1d0/(S*Color)                                                      
         !only scalars incoming=> no spin-averaging

! no flux and averaging factor!

      call QQ2xx_M2AmpInit2to2(1,channel,corrProp,corrVert,corrBox,QQ2xx_M2Tree,QQ2xx_M2Prop,QQ2xx_M2Vert,QQ2xx_M2Box)

      f(1) = jacobi*intfac*flux*avgfac*QQ2xx_M2Tree 
      f(2) = jacobi*intfac*flux*avgfac*QQ2xx_M2Tree
      f(3) = 0d0
      f(4) = 0d0
      f(5) = 0d0
      f(6) = 0d0
      f(7) = 0d0
      f(8) = 0d0
       
	end

c*************************************************************
      subroutine QQ2xx_Bremsstrahl23(intres, error)

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"

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
      external QQ2xx_integrand4       !Mo: Neu von mir analog zu ntnt2bb_final.F
      external vegas
      external cuhre
 
      flags = 0
      mineval = 1000
      maxeval = 500000
      maxiter = 1000
      nstart = maxeval/MAXITER
      nincrease = maxeval/(2*MAXITER)

      RELACCURACY = 5D-3 !1D-8
      ABSACCURACY = 1D-16  !1D-20

!         call cuhre(ndim, ncomp, QQ2xx_integrand4,
!      &      dble(RELACCURACY), dble(ABSACCURACY),
!      &      flags, mineval, maxeval,
!      &      key,
!      &      nregions, neval, fail, intres, error, prob)

c	intcuh(1) = intres(1)
c	ercuh(1) = error(1)
c	failcuh = fail
c	probcuh(1) = prob(1)
        

        call vegas(ndim, ncomp, QQ2xx_integrand4,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      flags, MINEVAL, MAXEVAL,
     &      nstart, nincrease,
     &      neval, fail, intres, error, prob)
 
       
c integration factor already includes the trivial integration over phi

	intfac23 = sqrtS**2/(32d0*(2*pi)**4)                                           

c* flux 

      if (fluxflag.eq.1) then                                                        
c* in units of pb                                                                      
!         flux = hbar_c2/(4d0*(kappa(sqrtS**2,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)
c* in units of GeV^-2
      flux = 1d0/(4d0*(kappa(sqrtS**2,m1**2,m2**2)/(2d0*sqrtS))*sqrtS)  

      else

c* in units of cm^3/sec
c	 flux = 2.9979d-26*hbar_c2/sqrtS**2
c* in units of GeV^-2
      flux = 1d0/sqrtS**2

      endif

!Symmetryfactor

        if((xtype.ne.0).and.(xtype.ne.4).and.(xtype.ne.8).and.((xtype.eq.ytype).and.(ytype.eq.ztype))) then

         S=6d0                                                                               !3*2 for three identical Gauge-Bosons in the final state (e.g. three Gluon final state)

        else if(((xtype.eq.ytype).or.(xtype.eq.ztype).or.(ytype.eq.ztype)).and.(xtype.ne.0).and.(xtype.ne.4).and.(xtype.ne.8)) then

         S=2d0                                                                               !Necessary, because phase-space-integration includes the points (p1,q2) and (p2,q1) (1,2 are the particles , p and q
                                                                                             !two momenta) which doublecounts the same process if particle 1 and 2 are identical.
        elseif((f3type.eq.f4type).and.(f3type.ne.0).and.(f3gen.eq.f4gen).and.(total_chargeout.ne.0)) then

         S=2d0                                                                               !covers the case of quark-quark or antiquark-antiquark final state                                                                       

        else

         S=1d0

        endif

!Coloraveraging

        if(((sftype1.eq.3).or.(sftype1.eq.4)).and.((sftype2.eq.3).or.(sftype2.eq.4))) then
         
         if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then      !Weak fermionic final states need to be reduced by an overcounting factor NC, which is already included in the 2->3-amplitudes
          Color = NC*NC
         else
          Color = NC
         endif
      
        else

         if(((f3type.eq.1).or.(f3type.eq.2)).and.((f4type.eq.1).or.(f4type.eq.2))) then      !Weak fermionic final states need to be reduced by an overcounting factor NC, which is already included in the 2->3-amplitudes
          Color = NC
         else
          Color = 1d0
         endif
         
        endif

!final averaging-factor 

	avgfac = 1d0/(S*Color)                                                               !only scalars incoming=> no spin-averaging

            
	intres(1) = avgfac*flux*intfac23*intres(1)
	error(1) = avgfac*flux*intfac23*error(1)

	end

c*************************************************************
	subroutine QQ2xx_Integrand4(ndim, x, ncomp, f)

        implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
	integer ndim, ncomp
	double precision x(*), f(*)

	double precision eta,costh,x1max,x1min,x2max,x2min
	double precision cosxi, sinxi, sinth
	double precision QQ2xx_M2Tree23
	double precision test
	
	external QQ2xx_PhSpace

! ===== Integration boundaries and Jacobi determinant =====

! dipole configuration -> gluon integrated out               !Dipoles not yet implemented => use phase-space slicing
!       cutoff = 0d0
!       call QQ2xx_PhSpace(x,muu3,muu5,muu4,muu1,muu2)                   

! soft-gluon configuration -> no gluon integrated out
        call QQ2xx_PhSpace(x,muu3,muu4,muu5,muu1,muu2,cosxi)  

! ===== Scalar products - dipole configuration =====

!         x1 = xx1
!         x2 = xx3    
!         x3 = xx2  

! ===== Scalar products - soft-gluon configuration =====

      x1 = xx1
      x2 = xx2
      x3 = xx3
	
! ===== Single contributions to matrix element =====

      if((dabs(cosxi).lt.0.999999999999999d0).and.(x1.gt.(2d0*cutoff/sqrtS)).and.(x2.gt.(2d0*cutoff/sqrtS)).and.(x3.gt.(2d0*cutoff/sqrtS))) then   !cut out numerical errors leading to a |cosxi|>1 oder cosxi=NAN (tests x1 and x2) and cut on potential 
                                                                                                                                                !photon/Gluon-energies leading to on-shell-propagators (soft-photon/Gluon)-> see also NeuQ2qx-code and associated paper page 13

      f(1) = jacobi*QQ2xx_M2Tree23()

      else

      f(1) = 0d0

      endif

	end

! ================= General Phase-space element======================
! input are reduced masses of final state (3x) and initial state (2x)

      subroutine QQ2xx_PhSpace(x,mmr1,mmr2,mmr3,mmra,mmrb,cosxi)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

      double precision x(*)
      double precision mmr1,mmr2,mmr3,mmra,mmrb
      double precision eta, costh,p1k3
      double precision cosxi,sinxi,sinth
      double precision xx3max,xx3min,xx1max,xx1min
      double precision taubar,sigmabar,muplus,muminus


! ***** Integration boundaries and Jacobi determinant (from HitchHiker) *****
            
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

      sinxi = dsqrt(1-cosxi**2)
      sinth = dsqrt(1-costh**2)

! *******************************Scalar products *******************************

      if(cutoff.lt.zerocheck) then  

!Dipol subtraction

          p1k1=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
          p1k2=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)  !sQuark-sQuark-scalar product

      else

!Phase space slicing
      
          p1k1=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx1-dsqrt(xx1**2-4*mmr1**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*(costh*cosxi - sinth*sinxi*dsin(eta)))
          p1k3=sqrtS**2/4d0*((1d0-mmrb**2+mmra**2)*xx3-dsqrt(xx3**2-4*mmr3**2)*dsqrt(1d0-2d0*(mmra**2+mmrb**2)+(mmra**2-mmrb**2)**2)*costh)  !sQuark-Gluon-scalar product 

          p1k2=(sqrtS**2/2d0)*(1d0-mmrb**2+mmra**2)-p1k1-p1k3                                                                                !p1k2=p1p1+p1p2-p1k1-p1k3

      endif

!******************************************************************************

      end


#include "QQ2xx_Couplings.F" 
#include "QQ2xx_Kinematics.F"
#include "QQ2xx_M2Loop.F"
#include "QQ2xx_M2HHgen.F"
#include "QQ2xx_M2VHgen.F"
#include "QQ2xx_M2HVgen.F"
#include "QQ2xx_M2VVgen.F"
#include "QQ2xx_M2FFbargen.F"
#include "QQ2xx_M2Init2to2.F"
#include "QQ2xx_M2AmpInit2to2.F"
#include "QQ2xx_M2AmpCorr.F"
#include "QQ2xx_M2AmpHH.F"
#include "QQ2xx_M2AmpHV.F"
#include "QQ2xx_M2AmpVH.F"
#include "QQ2xx_M2AmpVV.F"
#include "QQ2xx_M2AmpFF.F"
#include "QQ2xx_M2AmpFFbar.F"
#include "QQ2xx_M2Init2to3.F"
#include "QQ2xx_M2Tree23.F"
#include "QQ2xx_M2Tree23FFbar.F"
#include "QQ2xx_M2Tree23HH.F"
#include "QQ2xx_M2Tree23VH.F"
#include "QQ2xx_M2Tree23HV.F"
#include "QQ2xx_M2Tree23VV.F"
#include "QQ2xx_M2BoxHH.F"
#include "QQ2xx_M2BoxVH.F"
#include "QQ2xx_M2BoxHV.F"
#include "QQ2xx_M2BoxVV.F"
#include "QQ2xx_SoftGluon.F"
#include "QQ2xx_SetCounterterms.F"
#include "QQ2xx_M2CoulombCorrections.F"
#include "QQ2xx_M2DigestAmp.F"
#include "QQ2xx_M2AmpDigest.F"
#include "QQ2xx_M2_onshell_subtraction.F"