	function ststQQ_M2VertexT()
	
	implicit none

#include "../util/DMNLO_LoopIntegrals.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_GenCouplings.h" 

	integer i,j,k,l,n,m,neu,neuc
	complex*16 g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3L,g3R,g4L,g4R, g2Lc, g2Rc, g1Lc, g1Rc, g3Rc, g3Lc, g4Lc, g4Rc
	double precision Ml0,Ml1
	complex*16 vcAL,vcAR,vcBL,vcBR,vcCL,vcCR,vcDL,vcDR
	double precision Krondelta


* functions & subroutines
	double precision colorf, colori
	complex*16 ststQQ_M2VertexT
	external myCget, myBget

	ststQQ_M2VertexT = 0d0


!Propagator

	do i=1,4 !hier gluino und neutralino unterscheiden
		DenGauT(i)=(tman - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
	enddo
	DenGauT(5)=(tman - MCha(1)**2) !+ cI*MCha(1)*wCha(1)
	DenGauT(6)=(tman - MCha(2)**2) !+ cI*MCha(2)*wCha(2)

	DenGlT = (tman-MGl**2)
    
c****************************************************
c************ Loop external squark (ifs1,itt1,igen1) & quark (ftt1,fgen1) with Gluino in t-channel ***********************
c****************************************************

      vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0
	  
	  
c-------------------------------Vertex1 ( Gluon , Fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
	  do i=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluR(i,k,l,m,n)
	  g0R = GluL(i,k,l,m,n)
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2 = VSfSf(1,isf1,i,itt1,k,igen1,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
! 	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)	
	
	  colori = -1/(2d0*Nc)

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-tman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu3*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
      
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-tman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu3*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
!      
!      
	  vcBL = vcBL + colori*g1*g2*g0R*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcCL = vcCL + colori*g1*g2*(g0R*MGl*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGl*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)
	  
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo

c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluL(j,l,k,n,m)
	  g0R = GluR(j,l,k,n,m)
	  g1L = GluR(j,l,ftt1,n,fgen1)
	  g1R = GluL(j,l,ftt1,n,fgen1)
	  g2L = GluR(isf1,itt1,k,igen1,m)
	  g2R = GluL(isf1,itt1,k,igen1,m)
    	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = -1/(2d0*Nc)
	
 	  vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGl*muu3*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGl*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGl*(xc2+xc0)*g1L*g2L*g0R))
     
	  vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGl*muu3*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGl*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGl*(xc2+xc0)*g1R*g2R*g0L))


	  vcBL = vcBL - colori* (muu3*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	  vcBR = vcBR - colori* (muu3*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	  vcCL = vcCL - colori* (MGl*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	  vcCR = vcCL - colori* (MGl*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	  vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	  vcDR = vcDR + colori* xc1*g1L*g2R*g0L
! 
		  enddo
	     enddo
	    enddo
       enddo
	  enddo

c-------------------------------Vertex3 ( Sfermion (j,l,n),Gluino, Gluon (M0,M1,M2))	
	  do j=1,2
	   do l=3,4
		do n=1,3

	
	  g0 = GS
	  g1L = GluR(j,l,ftt1,n,fgen1)
	  g1R = GluL(j,l,ftt1,n,fgen1)
	  g2 = VSfSf(1,isf1,j,itt1,l,igen1,n)
		
	  Ml0 = MSf(j,l,n)
	  Ml1 = MGl
	  Ml2 = 0d0
	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0
	
 	  vcAL = vcAL + colori*g0*g2*(g1L*(-xb0-Ml0**2*xc0+muu1**2*sqrtS**2*xc2+muu3**2*sqrtS**2*(xc0-xc2)+Ml1*MGl*(xc2-xc0))+
     -			g1R*(-muu3*sqrtS*MGl*(xc0+xc1-xc2)-Ml1*muu3*sqrtS*(xc1+xc2-xc0)))
     
	  vcAR = vcAR + colori*g0*g2*(g1R*(-xb0-Ml0**2*xc0-muu3**2*sqrtS**2*xc1+muu1**2*sqrtS**2*xc2+muu3**2*sqrtS**2*(xc0-xc2+xc1)+Ml1*MGl*(xc2-xc0))+
     -			g1L*(-muu3*sqrtS*MGl*(xc0+xc1-xc2)-Ml1*muu3*sqrtS*(xc1+xc2-xc0)))

	  vcBL = vcBL + colori* g0*g2*(g1L*Ml1*(xc2-xc0)-g1R*muu3*sqrtS*(xc0+xc1-xc2))
!      
	  vcBR = vcBR + colori* g0*g2*(g1R*Ml1*(xc2-xc0)-g1L*muu3*sqrtS*(xc0+xc1-xc2))
	
	
	  vcCL = vcCL - colori* g0*g2*g1R*(MGl*(xc0+xc1-xc2)+Ml1*(xc1+xc2-xc0))
	
	  vcCR = vcCL - colori* g0*g2*g1L*(MGl*(xc0+xc1-xc2)+Ml1*(xc1+xc2-xc0))
	
	
	  vcDL = vcDL + colori* g0*g2*g1R*(xc0+xc1-xc2)      
	
	  vcDR = vcDR + colori* g0*g2*g1R*(xc0+xc1-xc2)

	    enddo
	   enddo
      enddo
	
c-------------------------------Vertex4 ( Fermion (l,n), Gluon, Gluino (M0,M1,M2))	

	  do l=3,4
	   do n=1,3
	
	  g0 = GS
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2L = GluR(isf1,itt1,l,igen1,n)
      g2R = GluL(isf1,itt1,l,igen1,n)
		
	  Ml0 = Mf(l,n)
	  Ml1 = 0d0
	  Ml2 = MGl
	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = Nc/2d0
	
 	  vcAL = vcAL - colori*g0*g1*(g2L*(4d0*xb0+4d0*Ml0**2*xc0+2d0*(-tman+muu3**2*sqrtS**2+muu1**2*sqrtS**2)*xc2+4d0*muu3**2*sqrtS**2*xc1+2d0*Ml0*muu3*sqrtS*(xc1+xc2+xc0)-2d0*MGl*Ml2*xc2)+
     -			g2R*(4d0*Ml0*Ml2*xc0+2d0*Ml2*muu3*sqrtS*(xc1+xc2)-2d0*Ml0*MGl*xc2))
     
	  vcAR = vcAR - colori*g0*g1*(g2R*(4d0*xb0+4d0*Ml0**2*xc0+2d0*(-tman+muu3**2*sqrtS**2+muu1**2*sqrtS**2)*(xc1)+4d0*muu1**2*sqrtS**2*xc2+2d0*Ml1*muu3*sqrtS*(xc1+xc2+xc0)-2d0*MGl*Ml2*xc2)+
     -			g2L*(4d0*Ml0*Ml2*xc0+2d0*Ml2*muu3*sqrtS*(xc1+xc2)-2d0*Ml0*MGl*xc2))

	  vcBL = vcBL - colori* g0*g1*(2d0*g2L*Ml2*xc2-2d0*g2R*Ml0*(xc0+xc2))
!      
	  vcBR = vcBR - colori* g0*g1*(2d0*g2R*Ml2*xc2-2d0*g2L*Ml0*(xc0+xc2))
	
	
	  vcCL = vcCL - colori* g0*g1*(-2d0*g2R*Ml2*(xc1+xc2)+2d0*Ml0*g2L*(xc1+xc2+xc0))
	
	  vcCR = vcCL - colori* g0*g1*(-2d0*g2L*Ml2*(xc1+xc2)+2d0*Ml0*g2R*(xc1+xc2+xc0))
	
	
	  vcDL = vcDL - colori* 0d0      
	
	  vcDR = vcDR - colori* 0d0
	
	   enddo
	  enddo

  	  vcAL = vcAL/(4d0*Pi)**2
  	  vcAR = vcAR/(4d0*Pi)**2
  	  vcBL = vcBL/(4d0*Pi)**2
  	  vcBR = vcBR/(4d0*Pi)**2
  	  vcCL = vcCL/(4d0*Pi)**2
  	  vcCR = vcCR/(4d0*Pi)**2
  	  vcDL = vcDL/(4d0*Pi)**2
  	  vcDR = vcDR/(4d0*Pi)**2

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c -----------------------Amplitude square VC Gluino t-channel with gluino T-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)


      colorf=(Nc**2-1)/4d0

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1Rc*(g2L*g2Rc*muu4*muugl*
     -          (2*muu3*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g2Rc*muugl*
     -          (muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2Lc*g2R*muu4*
     -          (-2*muu3*tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) - 
     -            (muu1**2 - tred)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2L*g2Lc*(muu3*(muu4**2 + tred)*
     -             (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (-(tred*(1 + tred)) + muu1**2*(muu4**2 + tred))*
     -             (vcAR - 2*muugl*sqrtS*vcBR) + 
     -            muu2**2*(-(muu3*
     -                  (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -                    sqrtS*tred*vcBL)) + 
     -               muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -               (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR))))+
     -        g1Lc*(g2L*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu2**2*muu3*vcAR - muu3*muu4**2*vcAR - 
     -            muu3*tred*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL - 
     -            sqrtS*tred*vcBL + muu3**2*sqrtS*tred*vcBL + 
     -            muu4**2*sqrtS*tred*vcBL + 
     -            2*muu3*muugl*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)
     -           + g2L*g2Rc*muu4*
     -          (-(muu1**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -                 sqrtS*tred*vcBL)) + 
     -            muu3**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) + 
     -            tred*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            2*muu3*tred*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2Lc*g2R*muu4*muugl*
     -          (muu1**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2R*g2Rc*(-(muu3**2*muu4**2*vcAL) - tred*vcAL - 
     -            tred**2*vcAL - muu3*muu4**2*muugl*vcAR - 
     -            muu3*muugl*tred*vcAR + 
     -            2*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muugl*sqrtS*tred*vcBL + 
     -            2*muugl*sqrtS*tred**2*vcBL - 
     -            muu1**2*(muu2**2 - muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3*muu4**2*muugl**2*sqrtS*vcBR + 
     -            muu3*muu4**2*sqrtS*tred*vcBR + 
     -            muu3*muugl**2*sqrtS*tred*vcBR + 
     -            muu3*sqrtS*tred**2*vcBR + 
     -            muu2**2*(muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -               tred*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -               muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -                  sqrtS*tred*vcBR))))))/(DenGlT*dconjg(DenGlT)))
!     

! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino t-channel with gluino u-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  
      colorf=-1d0/2d0*CF

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1Lc*(g2L*g2Rc*muu4*
     -          (-(muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) + muu2**2*muu3*vcAR - muu3**3*vcAR - 
     -            muu3*muu4**2*vcAR + 
     -            muugl**2*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL - 
     -            sqrtS*tred*vcBL + 2*muu3**2*sqrtS*tred*vcBL + 
     -            muu4**2*sqrtS*tred*vcBL - sqrtS*tred**2*vcBL + 
     -            muu1**2*(muu3*vcAR + sqrtS*tred*vcBL) + 
     -            2*muu3*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -             muugl*sqrtS*vcBR) + 
     -         g2L*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu2**2*muu3*vcAR - muu3*muu4**2*vcAR - 
     -            muu3*tred*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL - 
     -            sqrtS*tred*vcBL + muu3**2*sqrtS*tred*vcBL + 
     -            muu4**2*sqrtS*tred*vcBL + 
     -            2*muu3*muugl*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)
     -           + g2Lc*g2R*muu4*muugl*
     -          (muu1**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2R*g2Rc*(muu1**2*(muu2**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) - 
     -            (-1 + muu2**2 + muu4**2 - tred)*tred*
     -             (vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(muu4**2 + tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) + 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR))) + 
     -      g1Rc*(g2L*g2Rc*muu4*muugl*
     -          (2*muu3*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2L*g2Lc*(muu3**3*
     -             (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) + 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*(muu4**2 + tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2*(muu2**2 - tred) + 
     -               tred*(1 - muu2**2 - muu4**2 + tred))*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g2Rc*muugl*
     -          (muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2Lc*g2R*muu4*
     -          (-(muu3**3*vcAL) - muu3*muu4**2*vcAL + muugl*vcAR - 
     -            2*muu3**2*muugl*vcAR - muu4**2*muugl*vcAR + 
     -            muugl*tred*vcAR + 2*muu3**3*muugl*sqrtS*vcBL + 
     -            2*muu3*muu4**2*muugl*sqrtS*vcBL + 
     -            muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muugl**2*sqrtS*vcBR + 
     -            2*muu3**2*muugl**2*sqrtS*vcBR + 
     -            muu4**2*muugl**2*sqrtS*vcBR - sqrtS*tred*vcBR + 
     -            2*muu3**2*sqrtS*tred*vcBR + 
     -            muu4**2*sqrtS*tred*vcBR - 
     -            muugl**2*sqrtS*tred*vcBR - sqrtS*tred**2*vcBR + 
     -            muu1**2*(-(muugl*vcAR) + 
     -               muu3*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -               muugl**2*sqrtS*vcBR + sqrtS*tred*vcBR)))))/
     -  (DenGlT*dconjg(DenGlU)))


      do n=1,6
! 
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino t-channel with Neutralino T-channel -----------
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  colorf = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))


	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g2L*g3Lc*g4Rc*muu4*
     -       (-(muu1**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -              sqrtS*tred*vcBL)) + 
     -         muu3**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) + 
     -         tred*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) - 
     -         2*muu3*tred*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g2R*g3Rc*g4Lc*muu4*
     -       (-2*muu3*tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -            sqrtS*tred*vcBR) - 
     -         (muu1**2 - tred)*
     -          (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -            sqrtS*tred*vcBR)) + 
     -      g2L*g3Rc*g4Lc*(muu3*(muu4**2 + tred)*
     -          (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) - 
     -         muu3**2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         (-(tred*(1 + tred)) + muu1**2*(muu4**2 + tred))*
     -          (vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu2**2*(-(muu3*
     -               (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -                 sqrtS*tred*vcBL)) + 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR))) + 
     -      g2R*g3Lc*g4Rc*(-(muu3**2*muu4**2*vcAL) - tred*vcAL - 
     -         tred**2*vcAL - muu3*muu4**2*muugl*vcAR - 
     -         muu3*muugl*tred*vcAR + 
     -         2*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -         2*muugl*sqrtS*tred*vcBL + 2*muugl*sqrtS*tred**2*vcBL - 
     -         muu1**2*(muu2**2 - muu4**2 - tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3*muu4**2*muugl**2*sqrtS*vcBR + 
     -         muu3*muu4**2*sqrtS*tred*vcBR + 
     -         muu3*muugl**2*sqrtS*tred*vcBR + 
     -         muu3*sqrtS*tred**2*vcBR + 
     -         muu2**2*(muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR))) + 
     -      (g2L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu2**2*muu3*vcAR - muu3*muu4**2*vcAR - 
     -            muu3*tred*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL - 
     -            sqrtS*tred*vcBL + muu3**2*sqrtS*tred*vcBL + 
     -            muu4**2*sqrtS*tred*vcBL + 
     -            2*muu3*muugl*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)
     -           + g2L*g3Rc*g4Rc*muu4*
     -          (2*muu3*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g3Lc*g4Lc*muu4*
     -          (muu1**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2R*g3Rc*g4Rc*
     -          (muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)))*muupt(n)))/
     -  (DenGlT*dconjg(DenGauT(n))))

! 
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino t-channel with Neutralino U-channel -----------
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     
      colorf = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g2L*g3Lc*g4Rc*muu4*
     -       (-(muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -            vcAL) + muu3*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAR + 
     -         muugl**2*sqrtS*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -         sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*tred*
     -          vcBL + 2*muu3*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2)*muugl*sqrtS*
     -          vcBR) + g2L*g3Rc*g4Lc*
     -       (muu3**3*(-(muugl*vcAL) + 
     -            sqrtS*(muugl**2 + tred)*vcBL) + 
     -         muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -          (-(muugl*vcAL) + sqrtS*(muugl**2 + tred)*vcBL) - 
     -         muu3**2*(muu4**2 + tred)*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         (muu1**2*muu2**2 + tred - 
     -            (muu1**2 + muu2**2 + muu4**2)*tred + tred**2)*
     -          (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g2R*g3Lc*g4Rc*(muu1**2*(muu2**2 - tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) - 
     -         (-1 + muu2**2 + muu4**2 - tred)*tred*
     -          (vcAL - 2*muugl*sqrtS*vcBL) - 
     -         muu3**2*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3**3*(-(muugl*vcAR) + 
     -            sqrtS*(muugl**2 + tred)*vcBR) + 
     -         muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -          (-(muugl*vcAR) + sqrtS*(muugl**2 + tred)*vcBR)) + 
     -      g2R*g3Rc*g4Lc*muu4*
     -       (-(muu3**3*vcAL) - muu3*muu4**2*vcAL + muugl*vcAR - 
     -         2*muu3**2*muugl*vcAR - muu4**2*muugl*vcAR + 
     -         muugl*tred*vcAR + 2*muu3**3*muugl*sqrtS*vcBL + 
     -         2*muu3*muu4**2*muugl*sqrtS*vcBL + 
     -         muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         sqrtS*(-1 + 2*muu3**2 + muu4**2 - tred)*
     -          (muugl**2 + tred)*vcBR + 
     -         muu1**2*(-(muugl*vcAR) + 
     -            muu3*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            sqrtS*(muugl**2 + tred)*vcBR)) + 
     -      (g2L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu2**2*muu3*vcAR - muu3*(muu4**2 + tred)*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL + 
     -            (-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            2*muu3*muugl*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)
     -           + g2R*g3Lc*g4Lc*muu4*
     -          (-2*muu3*muugl*vcAR + 
     -            muu1**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3*sqrtS*(muugl**2 + tred)*vcBR) + 
     -         g2L*g3Rc*g4Rc*muu4*
     -          (-2*muu3*muugl*vcAL + 
     -            2*muu3*sqrtS*(muugl**2 + tred)*vcBL - 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g3Rc*g4Rc*
     -          (muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + 
     -               sqrtS*(muugl**2 + tred)*vcBR) + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + sqrtS*(muugl**2 + tred)*vcBR)))*
     -       muupu(n)))/(DenGlT*dconjg(DenGauU(n))))
!      
       enddo
! 
! !      
! !      
! !      
! !  ! ! ****************************************************
! ! ! ! ! ************ lower Loop external squark (ifs2,itt2,igen2) & quark (ftt2,fgen2) with Gluino in t-channel ***********************
! ! ! ! ! ****************************************************
! !       
	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0

! ! 	
! ! ! c-------------------------------Vertex1 ( Gluon , Fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
	  do i=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluLt1(i,k,l,m,n)
	  g0R = GluRt1(i,k,l,m,n)
	  g1 = BR(1,ftt2,l,fgen2,n)
	  g2 = VSfSf(1,isf2,i,itt2,k,igen2,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
! 	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)	
	
	  colori = -1/(2d0*Nc)
	  
! 	  vcAL = vcAL + colori*g1*g2*(-(g0R*muugl*sqrtS*
!      -       (-(Ml1*(2*xc0 + xc2)) + muu4*sqrtS*(xc0 + 2*xc1 + xc2))) +
!      -      g0L*(xb0 + Ml0**2*xc0 + 
!      -       sqrtS*(Ml1*muu4*(2*xc0 + xc1 + xc2) + 
!      -          sqrtS*(muu4**2*xc1 - tred*(xc0 + 2*xc1 + xc2) + 
!      -             muu2**2*(xc0 + 2*xc1 + 3*xc2)))))
!      
! 	  vcAR = vcAR + colori* g1*g2*(-(g0L*muugl*sqrtS*
!      -       (-(Ml1*(2*xc0 + xc2)) + muu4*sqrtS*(xc0 + 2*xc1 + xc2))) +
!      -     g0R*(xb0 + Ml0**2*xc0 + 
!      -       sqrtS*(Ml1*muu4*(2*xc0 + xc1 + xc2) + 
!      -          sqrtS*(muu4**2*xc1 - tred*(xc0 + 2*xc1 + xc2) + 
!      -             muu2**2*(xc0 + 2*xc1 + 3*xc2))))) 

      vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-tman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu4*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
      
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-tman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu4*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))

      vcBL = vcBL + colori*g1*g2*g0R*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcCL = vcCL + colori*g1*g2*(g0R*MGl*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGl*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)
	  
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo

! ! ! ! 
! ! ! ! c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

      g0L = GluRt1(j,l,k,n,m)
	  g0R = GluLt1(j,l,k,n,m)
	  g1L = GluLt1(j,l,ftt2,n,fgen2)
	  g1R = GluRt1(j,l,ftt2,n,fgen2)
	  g2L = GluLt1(isf2,itt2,k,igen2,m)
	  g2R = GluRt1(isf2,itt2,k,igen2,m)
	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = -1/(2d0*Nc)
	
 	  vcAL = vcAL +colori*g0R*muugl*sqrtS*(g1R*g2R*muu4*sqrtS*xc1 + 
     -     g1L*g2L*(Ml0*xc0 + (Ml0 + Ml2)*xc2)) - 
     -  g0L*(g1R*g2R*muu4*sqrtS*
     -      (Ml2*(xc1 + xc2) + Ml0*(xc0 + xc1 + xc2)) + 
     -     g1L*g2L*(xb0 + Ml0*(Ml0 + Ml2)*xc0 + 
     -        sqrtS**2*(muu4**2*xc1 + muu2**2*xc2)))
     
	  vcAR = vcAR +colori*(-(g0L*muugl*sqrtS*(g1L*g2L*muu4*sqrtS*xc1 + 
     -       g1R*g2R*(Ml0*xc0 + (Ml0 + Ml2)*xc2))) + 
     -  g0R*(g1L*g2L*muu4*sqrtS*
     -      (Ml2*(xc1 + xc2) + Ml0*(xc0 + xc1 + xc2)) + 
     -     g1R*g2R*(xb0 + Ml0*(Ml0 + Ml2)*xc0 + 
     -        sqrtS**2*(muu4**2*xc1 + muu2**2*xc2))))


	  vcBL = vcBL + colori*g0R*(g1R*g2R*muu4*sqrtS*xc1 + g1L*g2L*(Ml2*xc2 + Ml0*(xc0 + xc2)))
	  vcBR = vcBR + colori*(-(g0L*(g1L*g2L*muu4*sqrtS*xc1 + g1R*g2R*(Ml2*xc2 + Ml0*(xc0 + xc2)))))
	
	  vcCL = vcCL + colori*(-(g1R*g2R*(-(g0R*muugl*sqrtS*xc1) + g0L*Ml2*(xc1 + xc2) + 
     -      g0L*Ml0*(xc0 + xc1 + xc2))))
     
	  vcCR = vcCL + colori*g1L*g2L*(-(g0L*muugl*sqrtS*xc1) + g0R*Ml2*(xc1 + xc2) + 
     -    g0R*Ml0*(xc0 + xc1 + xc2))
	  
	  vcDL = vcDL + colori*g0R*g1R*g2R*xc1
	  vcDR = vcDR + colori*(-(g0L*g1L*g2L*xc1)) 
	  
		  enddo
	     enddo
	    enddo
       enddo
	  enddo

! c-------------------------------Vertex3 ( Sfermion (j,l,n),Gluino, Gluon (M0,M1,M2))	
	  do j=1,2
	   do l=3,4
		do n=1,3

	
	  g0 = GS
	  g1L = GluLt1(j,l,ftt2,n,fgen2)
	  g1R = GluRt1(j,l,ftt2,n,fgen2)
	  g2 = VSfSf(1,isf2,j,itt2,l,igen2,n)
		
	  Ml0 = MSf(j,l,n)
	  Ml1 = MGl
	  Ml2 = 0d0
	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0
	
 	  vcAL = vcAL + colori*g0*g2*(g1R*muu4*sqrtS*
     -     (muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2)) - 
     -    g1L*(xb0 + Ml0**2*xc0 - 
     -       sqrtS*(muu4**2*sqrtS*(xc0 - xc2) + muu2**2*sqrtS*xc2 + 
     -          Ml1*muugl*(-xc0 + xc2))))
 	  
	  vcAR = vcAR + colori*g0*g2*(g1L*muu4*sqrtS*
     -     (muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2)) - 
     -    g1R*(xb0 + Ml0**2*xc0 - 
     -       sqrtS*(muu4**2*sqrtS*(xc0 - xc2) + muu2**2*sqrtS*xc2 + 
     -          Ml1*muugl*(-xc0 + xc2))))
	  
	  vcBL = vcBL + colori*g0*g2*(g1R*muu4*sqrtS*(xc0 + xc1 - xc2) + g1L*Ml1*(-xc0 + xc2))
     
	  vcBR = vcBR + colori*g0*g2*(g1L*muu4*sqrtS*(xc0 + xc1 - xc2) + g1R*Ml1*(-xc0 + xc2))
	
	
	  vcCL = vcCL + colori*g0*g1R*g2*(muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2))
	  vcCR = vcCL + colori*g0*g1L*g2*(muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2))
	
	
	  vcDL = vcDL + colori*g0*g1R*g2*(xc0 + xc1 - xc2)
	  
	  vcDR = vcDR + colori*g0*g1L*g2*(xc0 + xc1 - xc2)

	    enddo
	   enddo
	  enddo


c-------------------------------Vertex4 ( Fermion (l,n), Gluon, Gluino (M0,M1,M2))	
! 
	  do l=3,4
	   do n=1,3
	
	  g0 = GS
	  g1 = BR(1,ftt2,l,fgen2,n)
	  g2L = GluLt1(isf2,itt2,l,igen2,n)
      g2R = GluRt1(isf2,itt2,l,igen2,n)
		
	  Ml0 = Mf(l,n)
	  Ml1 = 0d0
	  Ml2 = MGl
	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0

     
 	  vcAL = vcAL + colori*2*g0*g1*(g2R*Ml2*muu4*sqrtS*(xc1 + xc2) - 
     -    g2R*Ml0*(2*Ml2*xc0 + muugl*sqrtS*(xc0 + xc2)) + 
     -    g2L*(2*xb0 + 2*Ml0**2*xc0 - 
     -       Ml0*muu4*sqrtS*(xc0 + xc1 + xc2) + 
     -       sqrtS*(sqrtS*(muu2**2 + muu4**2 - tred)*xc1 + 
     -          Ml2*muugl*xc2 + 2*muu2**2*sqrtS*xc2)))
	   	  
	  vcAR = vcAR + colori*2*g0*g1*(g2L*Ml2*muu4*sqrtS*(xc1 + xc2) - 
     -    g2L*Ml0*(2*Ml2*xc0 + muugl*sqrtS*(xc0 + xc2)) + 
     -    g2R*(2*xb0 + 2*Ml0**2*xc0 - 
     -       Ml0*muu4*sqrtS*(xc0 + xc1 + xc2) + 
     -       sqrtS*(sqrtS*(muu2**2 + muu4**2 - tred)*xc1 + 
     -          Ml2*muugl*xc2 + 2*muu2**2*sqrtS*xc2)))
!      
	  vcBL = vcBL + colori*(-2*g0*g1*(-(g2L*Ml2*xc2) + g2R*Ml0*(xc0 + xc2)))
	  vcBR = vcBR + colori*(-2*g0*g1*(-(g2R*Ml2*xc2) + g2L*Ml0*(xc0 + xc2)))
	
	  vcCL = vcCL + colori*(-2*g0*g1*(-(g2R*Ml2*(xc1 + xc2)) + g2L*Ml0*(xc0 + xc1 + xc2)))
	  vcCR = vcCL + colori*(-2*g0*g1*(-(g2L*Ml2*(xc1 + xc2)) + g2R*Ml0*(xc0 + xc1 + xc2)))
	
	  vcDL = vcDL + colori* 0d0      
	
	  vcDR = vcDR + colori* 0d0
	
	   enddo
	  enddo

  	  vcAL = vcAL/(4d0*Pi)**2
  	  vcAR = vcAR/(4d0*Pi)**2
  	  vcBL = vcBL/(4d0*Pi)**2
  	  vcBR = vcBR/(4d0*Pi)**2
  	  vcCL = vcCL/(4d0*Pi)**2
  	  vcCR = vcCR/(4d0*Pi)**2
  	  vcDL = vcDL/(4d0*Pi)**2
  	  vcDR = vcDR/(4d0*Pi)**2

! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! c -----------------------Amplitude square VC Gluino t-channel with gluino T-channel -----------
! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)


      colorf= Nc/2d0*CF

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1R*(g1Lc*g2Lc*muu3*muugl*
     -          (-2*muu4*muugl*vcAR + 
     -            muu2**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBR) + 
     -         g1Lc*g2Rc*muu3*
     -          (-2*muu4*tred*vcAL + muu2**2*muugl*vcAR - 
     -            muu4**2*muugl*vcAR - muugl*tred*vcAR + 
     -            4*muu4*muugl*sqrtS*tred*vcBL + 
     -            sqrtS*(-muu2**2 + muu4**2 + tred)*(muugl**2 + tred)*
     -             vcBR) + g1Rc*g2Rc*muugl*
     -          (muugl*vcAR + 
     -            muu1**2*muu4*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4*(tred*vcAL + muu4*muugl*vcAR - 
     -               2*muugl*sqrtS*tred*vcBL) + 
     -            (-1 + muu4**2)*sqrtS*(muugl**2 + tred)*vcBR + 
     -            muu3**2*(-(muu4*vcAL) - muugl*vcAR + 
     -               2*muu4*muugl*sqrtS*vcBL + 
     -               sqrtS*(muugl**2 + tred)*vcBR)) + 
     -         g1Rc*g2Lc*(-(muu3**2*muu4**2*vcAL) - tred*vcAL - 
     -            tred**2*vcAL - muu3**2*muu4*muugl*vcAR - 
     -            muu4*muugl*tred*vcAR + 
     -            2*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muugl*sqrtS*tred*vcBL + 
     -            2*muugl*sqrtS*tred**2*vcBL + 
     -            muu2**2*(muu3**2 + tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu4*sqrtS*(muu3**2 + tred)*(muugl**2 + tred)*
     -             vcBR + muu1**2*
     -             (muu4*muugl*vcAR - 
     -               muu2**2*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -               muu4**2*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -               tred*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -               muu4*sqrtS*(muugl**2 + tred)*vcBR))) + 
     -      g1L*(g1Lc*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu1**2*muu4*vcAR - muu4*(muu3**2 + tred)*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL + 
     -            (-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            2*muu4*muugl*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)
     -           + g1Rc*g2Rc*muu3*muugl*
     -          (-2*muu4*muugl*vcAL + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBL - 
     -            muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1Rc*g2Lc*muu3*
     -          (-(muugl*tred*vcAL) + 
     -            sqrtS*tred*(muugl**2 + tred)*vcBL + 
     -            muu2**2*(muugl*vcAL - 
     -               sqrtS*(muugl**2 + tred)*vcBL) + 
     -            muu4**2*(-(muugl*vcAL) + 
     -               sqrtS*(muugl**2 + tred)*vcBL) - 
     -            2*muu4*tred*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1Lc*g2Rc*(muu3**2*
     -             (muu2**2*vcAR - muu4*(muugl*vcAL + muu4*vcAR) + 
     -               muu4*sqrtS*(muugl**2 + tred)*vcBL + 
     -               2*(-muu2 + muu4)*(muu2 + muu4)*muugl*sqrtS*vcBR)+
     -              tred*(-(muu4*muugl*vcAL) + 
     -               muu4*sqrtS*(muugl**2 + tred)*vcBL + 
     -               (-1 + muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)
     -               ) + muu1**2*
     -             (muu4*muugl*vcAL - 
     -               muu4*sqrtS*(muugl**2 + tred)*vcBL + 
     -               muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -               (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR))))))/
     -   (DenGlT*dconjg(DenGlT)))
! ! ! ! 	
! ! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! c -----------------------Amplitude square VC Gluino t-channel with gluino U-channel -----------
! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  
      colorf = -1d0/2d0*CF

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1L*(g1Rc*g2Lc*muu3*
     -          (-(muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) + 
     -            muu4*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAR + 
     -            muugl**2*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             tred*vcBL + 
     -            2*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -             muugl*sqrtS*vcBR) + 
     -         g1Lc*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu1**2*muu4*vcAR - muu4*(muu3**2 + tred)*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL + 
     -            (-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            2*muu4*muugl*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)
     -           + g1Rc*g2Rc*muu3*muugl*
     -          (-2*muu4*muugl*vcAL + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBL - 
     -            muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1Lc*g2Rc*(muu4**3*
     -             (-(muugl*vcAL) + sqrtS*(muugl**2 + tred)*vcBL) + 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (-(muugl*vcAL) + sqrtS*(muugl**2 + tred)*vcBL) - 
     -            muu4**2*(muu3**2 + tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2*muu2**2 + tred - 
     -               (muu1**2 + muu2**2 + muu3**2)*tred + tred**2)*
     -             (vcAR - 2*muugl*sqrtS*vcBR))) + 
     -      g1R*(g1Lc*g2Lc*muu3*muugl*
     -          (-2*muu4*muugl*vcAR + 
     -            muu2**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBR) + 
     -         g1Lc*g2Rc*muu3*
     -          (muu2**2*muu4*vcAL - muu3**2*muu4*vcAL - 
     -            muu4**3*vcAL + muugl*vcAR - muu2**2*muugl*vcAR - 
     -            muu3**2*muugl*vcAR - 2*muu4**2*muugl*vcAR + 
     -            muugl*tred*vcAR - 2*muu2**2*muu4*muugl*sqrtS*vcBL + 
     -            2*muu3**2*muu4*muugl*sqrtS*vcBL + 
     -            2*muu4**3*muugl*sqrtS*vcBL + 
     -            muu1**2*muu4*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             (muugl**2 + tred)*vcBR) + 
     -         g1Rc*g2Rc*muugl*
     -          (muugl*vcAR + 
     -            muu1**2*muu4*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4*(tred*vcAL + muu4*muugl*vcAR - 
     -               2*muugl*sqrtS*tred*vcBL) + 
     -            (-1 + muu4**2)*sqrtS*(muugl**2 + tred)*vcBR + 
     -            muu3**2*(-(muu4*vcAL) - muugl*vcAR + 
     -               2*muu4*muugl*sqrtS*vcBL + 
     -               sqrtS*(muugl**2 + tred)*vcBR)) + 
     -         g1Rc*g2Lc*(tred*vcAL - muu2**2*tred*vcAL - 
     -            muu4**2*tred*vcAL + tred**2*vcAL + 
     -            muu4*muugl*vcAR - muu4**3*muugl*vcAR + 
     -            muu4*muugl*tred*vcAR - 2*muugl*sqrtS*tred*vcBL + 
     -            2*muu2**2*muugl*sqrtS*tred*vcBL + 
     -            2*muu4**2*muugl*sqrtS*tred*vcBL - 
     -            2*muugl*sqrtS*tred**2*vcBL + 
     -            muu4*sqrtS*(-1 + muu4**2 - tred)*(muugl**2 + tred)*
     -             vcBR + muu1**2*
     -             (-(muu4*muugl*vcAR) + 
     -               (muu2**2 - tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -               muu4*sqrtS*(muugl**2 + tred)*vcBR) + 
     -            muu3**2*(-(muu4**2*vcAL) - tred*vcAL - 
     -               2*muu4*muugl*vcAR + 2*muu4**2*muugl*sqrtS*vcBL + 
     -               2*muugl*sqrtS*tred*vcBL + 
     -               2*muu4*sqrtS*(muugl**2 + tred)*vcBR)))))/
     -  (DenGlT*dconjg(DenGlU)))

        do n=1,6
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino t-channel with Neutralino T-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  colorf = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))


	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g2L*g3Lc*g4Rc*muu4*
     -       (-(muu1**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -              sqrtS*tred*vcBL)) + 
     -         muu3**2*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) + 
     -         tred*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) - 
     -         2*muu3*tred*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g2R*g3Rc*g4Lc*muu4*
     -       (-2*muu3*tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -            sqrtS*tred*vcBR) - 
     -         (muu1**2 - tred)*
     -          (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -            sqrtS*tred*vcBR)) + 
     -      g2L*g3Rc*g4Lc*(muu3*(muu4**2 + tred)*
     -          (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -            sqrtS*tred*vcBL) - 
     -         muu3**2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         (-(tred*(1 + tred)) + muu1**2*(muu4**2 + tred))*
     -          (vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu2**2*(-(muu3*
     -               (-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -                 sqrtS*tred*vcBL)) + 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR))) + 
     -      g2R*g3Lc*g4Rc*(-(muu3**2*muu4**2*vcAL) - tred*vcAL - 
     -         tred**2*vcAL - muu3*muu4**2*muugl*vcAR - 
     -         muu3*muugl*tred*vcAR + 
     -         2*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -         2*muugl*sqrtS*tred*vcBL + 2*muugl*sqrtS*tred**2*vcBL - 
     -         muu1**2*(muu2**2 - muu4**2 - tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3*muu4**2*muugl**2*sqrtS*vcBR + 
     -         muu3*muu4**2*sqrtS*tred*vcBR + 
     -         muu3*muugl**2*sqrtS*tred*vcBR + 
     -         muu3*sqrtS*tred**2*vcBR + 
     -         muu2**2*(muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR))) + 
     -      (g2L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu2**2*muu3*vcAR - muu3*muu4**2*vcAR - 
     -            muu3*tred*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL - 
     -            sqrtS*tred*vcBL + muu3**2*sqrtS*tred*vcBL + 
     -            muu4**2*sqrtS*tred*vcBL + 
     -            2*muu3*muugl*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)
     -           + g2L*g3Rc*g4Rc*muu4*
     -          (2*muu3*(-(muugl*vcAL) + muugl**2*sqrtS*vcBL + 
     -               sqrtS*tred*vcBL) - 
     -            muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g3Lc*g4Lc*muu4*
     -          (muu1**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)) + 
     -         g2R*g3Rc*g4Rc*
     -          (muu2**2*muu3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**2*(-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR) + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + muugl**2*sqrtS*vcBR + 
     -               sqrtS*tred*vcBR)))*muupt(n)))/
     -  (DenGlT*dconjg(DenGauT(n))))

! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino t-channel with Neutralino U-channel -----------
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      
! 	  
	  colorf = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1L*g3Rc*g4Lc*muu3*
     -       (-(muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -            vcAL) + muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAR + 
     -         muugl**2*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -         sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*
     -          vcBL + 2*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2)*muugl*sqrtS*
     -          vcBR) + g1R*g3Lc*g4Rc*muu3*
     -       (muu2**2*muu4*vcAL - muu3**2*muu4*vcAL - muu4**3*vcAL + 
     -         muugl*vcAR - muu2**2*muugl*vcAR - muu3**2*muugl*vcAR - 
     -         2*muu4**2*muugl*vcAR + muugl*tred*vcAR - 
     -         2*muu2**2*muu4*muugl*sqrtS*vcBL + 
     -         2*muu3**2*muu4*muugl*sqrtS*vcBL + 
     -         2*muu4**3*muugl*sqrtS*vcBL + 
     -         muu1**2*muu4*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -          (muugl**2 + tred)*vcBR) + 
     -      g1L*g3Lc*g4Rc*(muu4**3*
     -          (-(muugl*vcAL) + sqrtS*(muugl**2 + tred)*vcBL) + 
     -         muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -          (-(muugl*vcAL) + sqrtS*(muugl**2 + tred)*vcBL) - 
     -         muu4**2*(muu3**2 + tred)*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         (muu1**2*muu2**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2)*tred + tred**2)*
     -          (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g1R*g3Rc*g4Lc*(tred*vcAL - muu2**2*tred*vcAL - 
     -         muu4**2*tred*vcAL + tred**2*vcAL + muu4*muugl*vcAR - 
     -         muu4**3*muugl*vcAR + muu4*muugl*tred*vcAR - 
     -         2*muugl*sqrtS*tred*vcBL + 
     -         2*muu2**2*muugl*sqrtS*tred*vcBL + 
     -         2*muu4**2*muugl*sqrtS*tred*vcBL - 
     -         2*muugl*sqrtS*tred**2*vcBL + 
     -         muu4*sqrtS*(-1 + muu4**2 - tred)*(muugl**2 + tred)*
     -          vcBR + muu1**2*
     -          (-(muu4*muugl*vcAR) + 
     -            (muu2**2 - tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu4*sqrtS*(muugl**2 + tred)*vcBR) + 
     -         muu3**2*(-(muu4**2*vcAL) - tred*vcAL - 
     -            2*muu4*muugl*vcAR + 2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muugl*sqrtS*tred*vcBL + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBR)) + 
     -      (g1L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) + 
     -            muu1**2*muu4*vcAR - muu4*(muu3**2 + tred)*vcAR + 
     -            (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*vcBL + 
     -            (-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            2*muu4*muugl*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)
     -           + g1R*g3Lc*g4Lc*muu3*
     -          (-2*muu4*muugl*vcAR + 
     -            muu2**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4**2*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            tred*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBR) + 
     -         g1L*g3Rc*g4Rc*muu3*
     -          (-2*muu4*muugl*vcAL + 
     -            2*muu4*sqrtS*(muugl**2 + tred)*vcBL - 
     -            muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1R*g3Rc*g4Rc*
     -          (muugl*vcAR + 
     -            muu1**2*muu4*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4*(tred*vcAL + muu4*muugl*vcAR - 
     -               2*muugl*sqrtS*tred*vcBL) + 
     -            (-1 + muu4**2)*sqrtS*(muugl**2 + tred)*vcBR + 
     -            muu3**2*(-(muu4*vcAL) - muugl*vcAR + 
     -               2*muu4*muugl*sqrtS*vcBL + 
     -               sqrtS*(muugl**2 + tred)*vcBR)))*muupu(n)))/
     -  (DenGlT*dconjg(DenGauU(n))))     



        enddo

        
        
        
! ! ! ! ! ************************************************************************************************************************************************************
! ! ! ! ! ************ Upper Loop external squark (ifs2,itt2,igen2) & quark (ftt1,fgen1) with Gluino in u-channel ****************************************************
! ! ! ! ! ************************************************************************************************************************************************************
! ! ! ! !       
	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0

c-------------------------------Vertex1 ( Gluon , Fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))---------------------------------------------------------------
	  do i=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluR(i,k,l,m,n)
	  g0R = GluL(i,k,l,m,n)
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2 = VSfSf(1,isf2,i,itt2,k,igen2,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
! 	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)	
	
	  colori = -1/(2d0*Nc)

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-uman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu3*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
      
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-uman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu3*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
!      
!      
	  vcBL = vcBL + colori*g1*g2*g0R*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcCL = vcCL + colori*g1*g2*(g0R*MGl*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGl*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)
	  
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo


! ! ! 
! c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluL(j,l,k,n,m)
	  g0R = GluR(j,l,k,n,m)
	  g1L = GluR(j,l,ftt1,n,fgen1)
	  g1R = GluL(j,l,ftt1,n,fgen1)
	  g2L = GluR(isf2,itt2,k,igen2,m)
	  g2R = GluL(isf2,itt2,k,igen2,m)
    	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = -1/(2d0*Nc)
	
 	  vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGl*muu3*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGl*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGl*(xc2+xc0)*g1L*g2L*g0R))
     
	  vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGl*muu3*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGl*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGl*(xc2+xc0)*g1R*g2R*g0L))


	  vcBL = vcBL - colori* (muu3*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	  vcBR = vcBR - colori* (muu3*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	  vcCL = vcCL - colori* (MGl*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	  vcCR = vcCL - colori* (MGl*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	  vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	  vcDR = vcDR + colori* xc1*g1L*g2R*g0L
! 
		  enddo
	     enddo
	    enddo
       enddo
	  enddo

	  

c-------------------------------Vertex3 ( Sfermion (j,l,n),Gluino, Gluon (M0,M1,M2))	
	  do j=1,2
	   do l=3,4
		do n=1,3

	
	  g0 = GS
	  g1L = GluR(j,l,ftt1,n,fgen1)
	  g1R = GluL(j,l,ftt1,n,fgen1)
	  g2 = VSfSf(1,isf2,j,itt2,l,igen2,n)
		
	  Ml0 = MSf(j,l,n)
	  Ml1 = MGl
	  Ml2 = 0d0
	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0
	
 	  vcAL = vcAL + colori*g0*g2*(g1L*(-xb0-Ml0**2*xc0+muu2**2*sqrtS**2*xc2+muu3**2*sqrtS**2*(xc0-xc2)+Ml1*MGl*(xc2-xc0))+
     -			g1R*(-muu3*sqrtS*MGl*(xc0+xc1-xc2)-Ml1*muu3*sqrtS*(xc1+xc2-xc0)))
     
	  vcAR = vcAR + colori*g0*g2*(g1R*(-xb0-Ml0**2*xc0-muu3**2*sqrtS**2*xc1+muu2**2*sqrtS**2*xc2+muu3**2*sqrtS**2*(xc0-xc2+xc1)+Ml1*MGl*(xc2-xc0))+
     -			g1L*(-muu3*sqrtS*MGl*(xc0+xc1-xc2)-Ml1*muu3*sqrtS*(xc1+xc2-xc0)))

	  vcBL = vcBL + colori* g0*g2*(g1L*Ml1*(xc2-xc0)-g1R*muu3*sqrtS*(xc0+xc1-xc2))
!      
	  vcBR = vcBR + colori* g0*g2*(g1R*Ml1*(xc2-xc0)-g1L*muu3*sqrtS*(xc0+xc1-xc2))
	
	
	  vcCL = vcCL - colori* g0*g2*g1R*(MGl*(xc0+xc1-xc2)+Ml1*(xc1+xc2-xc0))
	
	  vcCR = vcCL - colori* g0*g2*g1L*(MGl*(xc0+xc1-xc2)+Ml1*(xc1+xc2-xc0))
	
	
	  vcDL = vcDL + colori* g0*g2*g1R*(xc0+xc1-xc2)      
	
	  vcDR = vcDR + colori* g0*g2*g1R*(xc0+xc1-xc2)

	    enddo
	   enddo
	  enddo

	
c-------------------------------Vertex4 ( Fermion (l,n), Gluon, Gluino (M0,M1,M2))	

	  do l=3,4
	   do n=1,3
	
	  g0 = GS
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2L = GluR(isf2,itt2,l,igen2,n)
      g2R = GluL(isf2,itt2,l,igen2,n)
		
	  Ml0 = Mf(l,n)
	  Ml1 = 0d0
	  Ml2 = MGl
	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = Nc/2d0
	
 	  vcAL = vcAL - colori*g0*g1*(g2L*(4d0*xb0+4d0*Ml0**2*xc0+2d0*(-uman+muu3**2*sqrtS**2+muu2**2*sqrtS**2)*
     -          xc2+4d0*muu3**2*sqrtS**2*xc1+2d0*Ml0*muu3*sqrtS*(xc1+xc2+xc0)-2d0*MGl*Ml2*xc2)+
     -			g2R*(4d0*Ml0*Ml2*xc0+2d0*Ml2*muu3*sqrtS*(xc1+xc2)-2d0*Ml0*MGl*xc2))
     
	  vcAR = vcAR - colori*g0*g1*(g2R*(4d0*xb0+4d0*Ml0**2*xc0+2d0*(-uman+muu3**2*sqrtS**2+muu2**2*sqrtS**2)*
     -          (xc1)+4d0*muu2**2*sqrtS**2*xc2+2d0*Ml1*muu3*sqrtS*(xc1+xc2+xc0)-2d0*MGl*Ml2*xc2)+
     -			g2L*(4d0*Ml0*Ml2*xc0+2d0*Ml2*muu3*sqrtS*(xc1+xc2)-2d0*Ml0*MGl*xc2))

	  vcBL = vcBL - colori* g0*g1*(2d0*g2L*Ml2*xc2-2d0*g2R*Ml0*(xc0+xc2))
!      
	  vcBR = vcBR - colori* g0*g1*(2d0*g2R*Ml2*xc2-2d0*g2L*Ml0*(xc0+xc2))
	
	
	  vcCL = vcCL - colori* g0*g1*(-2d0*g2R*Ml2*(xc1+xc2)+2d0*Ml0*g2L*(xc1+xc2+xc0))
	
	  vcCR = vcCL - colori* g0*g1*(-2d0*g2L*Ml2*(xc1+xc2)+2d0*Ml0*g2R*(xc1+xc2+xc0))
	
	
	  vcDL = vcDL - colori* 0d0      
	
	  vcDR = vcDR - colori* 0d0
	
	   enddo
	  enddo

  	  vcAL = vcAL/(4d0*Pi)**2
  	  vcAR = vcAR/(4d0*Pi)**2
  	  vcBL = vcBL/(4d0*Pi)**2
  	  vcBR = vcBR/(4d0*Pi)**2
  	  vcCL = vcCL/(4d0*Pi)**2
  	  vcCR = vcCR/(4d0*Pi)**2
  	  vcDL = vcDL/(4d0*Pi)**2
  	  vcDR = vcDR/(4d0*Pi)**2
  	  
! ! ! ! ! 
! ! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! c -----------------------Amplitude square VC Gluino u-channel with gluino u-channel -----------
! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  colorf = Nc/2d0*CF
    
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1Lc*(g2L*g2Rc*muu4*
     -          (-(muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) + 
     -            muugl**2*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -             (-2*muu3*vcAR + 
     -               sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            4*muu3*muugl*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g2L*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu3*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) +
     -           g2Lc*g2R*muu4*muugl*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) - 2*muu3*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g2R*g2Rc*(-(muu1**2*muu2**2*vcAL) - 
     -            muu1**2*muu4**2*vcAL + muu1**2*tred*vcAL - 
     -            muu3**3*muugl*vcAR + 
     -            muu3*muugl*(1 - muu2**2 - 2*muu4**2 + tred)*vcAR + 
     -            2*muu1**2*muu2**2*muugl*sqrtS*vcBL + 
     -            2*muu1**2*muu4**2*muugl*sqrtS*vcBL - 
     -            2*muu1**2*muugl*sqrtS*tred*vcBL - 
     -            muu3**4*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) - 
     -            (muu4**4 + tred*(1 - muu2**2 + tred) - 
     -               muu4**2*(1 + 2*tred))*(vcAL - 2*muugl*sqrtS*vcBL)
     -              + muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR)) + 
     -      g1Rc*(g2L*g2Rc*muu4*muugl*
     -          (-2*muu3*muugl*vcAL + 2*muu3**3*sqrtS*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu4**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu1**2 + muu4**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2L*g2Lc*(muu3**5*sqrtS*vcBL + 
     -            muu3**3*(-(muugl*vcAL) + 
     -               sqrtS*(-2 + muu1**2 + 2*muu2**2 + 3*muu4**2 + 
     -                  muugl**2 - 2*tred)*vcBL) + 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (-(muugl*vcAL) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBL) - 
     -            muu3**4*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-muu4**2 + (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2Lc*g2R*muu4*
     -          (-2*muu3**3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            2*muu3*(-1 + muu1**2 + muu2**2 + muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu3**4*sqrtS*vcBR + 
     -            muu3**2*(-2*muugl*vcAR + 
     -               sqrtS*(-3 + 3*muu1**2 + 2*muu2**2 + 3*muu4**2 + 
     -                  2*muugl**2 - 3*tred)*vcBR) + 
     -            (-1 + muu1**2 + muu4**2 - tred)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR)) + 
     -         g2R*g2Rc*muugl*
     -          (-(muu3**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu3**4*sqrtS*vcBR + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu3**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu4**2 + 
     -                  muugl**2 - tred)*vcBR)))))/
     -  (DenGlU*dconjg(DenGlU)))        
       
! !  ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! ! c -----------------------Amplitude square VC Gluino u-channel with gluino t-channel -----------
! ! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! 
	  colorf = -1d0/2d0*CF
	
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1Lc*(g2L*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu3*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR)+
     -          g2Lc*g2R*muu4*muugl*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) - 2*muu3*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g2R*g2Rc*(-((-1 + muu2**2 + muu4**2 - tred)*tred*
     -               (vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu3**2*(muu4**2 + tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3**3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR - 
     -            muu3*(muu2**2 - muu4**2 - tred)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu2**2 + muu4**2 + muugl**2 - tred)*
     -                vcBR) + 
     -            muu1**2*((muu2**2 - tred)*
     -                (vcAL - 2*muugl*sqrtS*vcBL) + 
     -               muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)) + 
     -         g2L*g2Rc*muu4*
     -          (-(muugl*tred*vcAL) - muu1**4*sqrtS*vcBL + 
     -            muu3**4*sqrtS*vcBL + 
     -            sqrtS*(-1 + muu2**2 + muu4**2 + muugl**2 - tred)*
     -             tred*vcBL + 
     -            muu3**2*(-(muugl*vcAL) + 
     -               (-1 + muu2**2 + muu4**2 + muugl**2)*sqrtS*vcBL) -
     -              muu3**3*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            muu3*(muu2 - muu4)*(muu2 + muu4)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) + 
     -            muu1**2*(muu3*vcAR - muugl**2*sqrtS*vcBL - 
     -               sqrtS*(-1 + muu2**2 + muu4**2 - 2*tred)*vcBL + 
     -               muugl*(vcAL - 2*muu3*sqrtS*vcBR)))) + 
     -      g1Rc*(g2Lc*g2R*muu4*
     -          (muu1**2*muu3*vcAL + muu2**2*muu3*vcAL - 
     -            muu3**3*vcAL - muu3*muu4**2*vcAL + 
     -            muu1**2*muugl*vcAR - muu3**2*muugl*vcAR - 
     -            muugl*tred*vcAR - 2*muu1**2*muu3*muugl*sqrtS*vcBL - 
     -            2*muu2**2*muu3*muugl*sqrtS*vcBL + 
     -            2*muu3**3*muugl*sqrtS*vcBL + 
     -            2*muu3*muu4**2*muugl*sqrtS*vcBL + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*(-muu1**2 + muu3**2 + tred)*vcBR
     -            ) + g2L*g2Rc*muu4*muugl*
     -          (-2*muu3*muugl*vcAL + 2*muu3**3*sqrtS*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu4**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu1**2 + muu4**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g2Rc*muugl*
     -          (-(muu3**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu3**4*sqrtS*vcBR + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu3**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu4**2 + 
     -                  muugl**2 - tred)*vcBR)) + 
     -         g2L*g2Lc*(-(muu2**4*muu3*sqrtS*vcBL) + 
     -            muu3**3*sqrtS*(muu4**2 + tred)*vcBL + 
     -            muu3*(muu4**2 + tred)*
     -             (-(muugl*vcAL) + 
     -               sqrtS*(-1 + muu1**2 + muu4**2 + muugl**2 - tred)*
     -                vcBL) + 
     -            tred*(1 - muu1**2 - muu4**2 + tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) - 
     -            muu3**2*(muu4**2 + tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR) + 
     -            muu2**2*(muu3*muugl*vcAL - muu3**3*sqrtS*vcBL - 
     -               muu3*sqrtS*(-1 + muu1**2 + muugl**2 - 2*tred)*
     -                vcBL + 
     -               (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR))))))
     -   /(DenGlU*dconjg(DenGlT)))        
!                 
!                 
         do n = 1,6       
 ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! c -----------------------Amplitude square VC Gluino u-channel with neutralino t-channel -----------
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !         
      colorf = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
                
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g2R*g3Rc*g4Lc*muu4*
     -       (muu1**2*muu3*vcAL + muu2**2*muu3*vcAL - muu3**3*vcAL - 
     -         muu3*muu4**2*vcAL + muu1**2*muugl*vcAR - 
     -         muu3**2*muugl*vcAR - muugl*tred*vcAR - 
     -         2*muu1**2*muu3*muugl*sqrtS*vcBL - 
     -         2*muu2**2*muu3*muugl*sqrtS*vcBL + 
     -         2*muu3**3*muugl*sqrtS*vcBL + 
     -         2*muu3*muu4**2*muugl*sqrtS*vcBL + 
     -         sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -            muugl**2 - tred)*(-muu1**2 + muu3**2 + tred)*vcBR) +
     -        g2R*g3Lc*g4Rc*
     -       (-((-1 + muu2**2 + muu4**2 - tred)*tred*
     -            (vcAL - 2*muugl*sqrtS*vcBL)) - 
     -         muu3**2*(muu4**2 + tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3**3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR - 
     -         muu3*(muu2**2 - muu4**2 - tred)*
     -          (-(muugl*vcAR) + 
     -            sqrtS*(-1 + muu2**2 + muu4**2 + muugl**2 - tred)*
     -             vcBR) + muu1**2*
     -          ((muu2**2 - tred)*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -            muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)) + 
     -      g2L*g3Lc*g4Rc*muu4*
     -       (-(muugl*tred*vcAL) - muu1**4*sqrtS*vcBL + 
     -         muu3**4*sqrtS*vcBL + 
     -         sqrtS*(-1 + muu2**2 + muu4**2 + muugl**2 - tred)*tred*
     -          vcBL + muu3**2*
     -          (-(muugl*vcAL) + 
     -            (-1 + muu2**2 + muu4**2 + muugl**2)*sqrtS*vcBL) - 
     -         muu3**3*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu3*(muu2 - muu4)*(muu2 + muu4)*
     -          (vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu1**2*(muu3*vcAR - muugl**2*sqrtS*vcBL - 
     -            sqrtS*(-1 + muu2**2 + muu4**2 - 2*tred)*vcBL + 
     -            muugl*(vcAL - 2*muu3*sqrtS*vcBR))) + 
     -      g2L*g3Rc*g4Lc*(-(muu2**4*muu3*sqrtS*vcBL) + 
     -         muu3**3*sqrtS*(muu4**2 + tred)*vcBL + 
     -         muu3*(muu4**2 + tred)*
     -          (-(muugl*vcAL) + 
     -            sqrtS*(-1 + muu1**2 + muu4**2 + muugl**2 - tred)*
     -             vcBL) + tred*(1 - muu1**2 - muu4**2 + tred)*
     -          (vcAR - 2*muugl*sqrtS*vcBR) - 
     -         muu3**2*(muu4**2 + tred)*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu2**2*(muu3*muugl*vcAL - muu3**3*sqrtS*vcBL - 
     -            muu3*sqrtS*(-1 + muu1**2 + muugl**2 - 2*tred)*
     -             vcBL + (muu1**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)
     -            )) + (g2L*g3Lc*g4Lc*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu3*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) +
     -           g2R*g3Lc*g4Lc*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) - 2*muu3*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g2L*g3Rc*g4Rc*muu4*
     -          (-2*muu3*muugl*vcAL + 2*muu3**3*sqrtS*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu4**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu1**2 + muu4**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g3Rc*g4Rc*
     -          (-(muu3**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu3**4*sqrtS*vcBR + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu3**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu4**2 + 
     -                  muugl**2 - tred)*vcBR)))*muupt(n)))/
     -  (DenGlU*dconjg(DenGauT(n))))                
!                 
!                 
!  ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! c -----------------------Amplitude square VC Gluino u-channel with neutralino u-channel -----------
! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !         
      colorf = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  
                
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g2L*g3Lc*g4Rc*muu4*
     -       (-(muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -            vcAL) + muugl**2*sqrtS*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -         (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-2*muu3*vcAR + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcBL) + 4*muu3*muugl*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g2L*g3Rc*g4Lc*
     -       (muu3**5*sqrtS*vcBL + 
     -         muu3**3*(-(muugl*vcAL) + 
     -            sqrtS*(-2 + muu1**2 + 2*muu2**2 + 3*muu4**2 + 
     -               muugl**2 - 2*tred)*vcBL) + 
     -         muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -          (-(muugl*vcAL) + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBL) - 
     -         muu3**4*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -         muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred)*
     -          (vcAR - 2*muugl*sqrtS*vcBR) - 
     -         (-muu4**2 + (muu4**2 - tred)**2 + 
     -            muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -            muu2**2*tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g2R*g3Lc*g4Rc*(-(muu3**4*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -         muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) - 
     -         (muu4**4 + tred*(1 - muu2**2 + tred) - 
     -            muu4**2*(1 + 2*tred))*(vcAL - 2*muugl*sqrtS*vcBL) + 
     -         muu3**5*sqrtS*vcBR + 
     -         muu3**3*(-(muugl*vcAR) + 
     -            sqrtS*(-2 + 2*muu2**2 + 3*muu4**2 + muugl**2 - 
     -               2*tred)*vcBR) + 
     -         muu1**2*(-((muu2**2 + muu4**2 - tred)*
     -               (vcAL - 2*muugl*sqrtS*vcBL)) + 
     -            muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR)+
     -           muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -          (-(muugl*vcAR) + 
     -            sqrtS*(-1 + muu2**2 + muu4**2 + muugl**2 - tred)*
     -             vcBR)) + g2R*g3Rc*g4Lc*muu4*
     -       (-2*muu3**3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -         2*muu3*(-1 + muu1**2 + muu2**2 + muu4**2 - tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) + 2*muu3**4*sqrtS*vcBR + 
     -         muu3**2*(-2*muugl*vcAR + 
     -            sqrtS*(-3 + 3*muu1**2 + 2*muu2**2 + 3*muu4**2 + 
     -               2*muugl**2 - 3*tred)*vcBR) + 
     -         (-1 + muu1**2 + muu4**2 - tred)*
     -          (-(muugl*vcAR) + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR)) + 
     -      (g2L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu3*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) +
     -          g2R*g3Lc*g4Lc*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) - 2*muu3*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g2L*g3Rc*g4Rc*muu4*
     -          (-2*muu3*muugl*vcAL + 2*muu3**3*sqrtS*vcBL + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu4**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu3**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu1**2 + muu4**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g2R*g3Rc*g4Rc*
     -          (-(muu3**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu3*(-1 + muu2**2 + 2*muu4**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu3**4*sqrtS*vcBR + 
     -            (-1 + muu4**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu3**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu4**2 + 
     -                  muugl**2 - tred)*vcBR)))*muupu(n)))/
     -  (DenGlU*dconjg(DenGauU(n))))        
                
                
         enddo       

         
         
         
         
         
         
         
         
         
         
         
c************************************************************************************************************************************************************
c************ Lower Loop external squark (ifs1,itt1,igen1) & quark (ftt2,fgen2) with Gluino in u-channel ****************************************************
c************************************************************************************************************************************************************
     
	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0

c--------------------------------------------------Vertex1 ( Gluon , Fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))-----------------------------------------------

	  do i=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

	  g0L = GluR(i,k,l,m,n)
	  g0R = GluL(i,k,l,m,n)
	  g1 = BR(1,ftt2,l,fgen2,n)
	  g2 = VSfSf(1,isf1,i,itt1,k,igen1,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
! 	
	  call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)	
	
	  colori = -1/(2d0*Nc)

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-uman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu4*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))
      
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-uman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu4*sqrtS*MGl*(2*xc1+xc2+2*xc0)-MGl*Ml1*(xc2+2*xc0)))

	  vcBL = vcBL + colori*g1*g2*g0R*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcCL = vcCL + colori*g1*g2*(g0R*MGl*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGl*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)
	  
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo

! ! ! c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3

      g0L = GluLt2(j,l,k,n,m)
	  g0R = GluRt2(j,l,k,n,m)
	  g1L = GluLt1(j,l,ftt1,n,fgen1)
	  g1R = GluRt1(j,l,ftt1,n,fgen1)
	  g2L = GluLt1(isf1,itt1,k,igen1,m)
	  g2R = GluRt1(isf1,itt1,k,igen1,m)
	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = -1/(2d0*Nc)
	
 	  vcAL = vcAL +colori*g0R*muugl*sqrtS*(g1R*g2R*muu4*sqrtS*xc1 + 
     -     g1L*g2L*(Ml0*xc0 + (Ml0 + Ml2)*xc2)) - 
     -  g0L*(g1R*g2R*muu4*sqrtS*
     -      (Ml2*(xc1 + xc2) + Ml0*(xc0 + xc1 + xc2)) + 
     -     g1L*g2L*(xb0 + Ml0*(Ml0 + Ml2)*xc0 + 
     -        sqrtS**2*(muu4**2*xc1 + muu1**2*xc2)))
     
	  vcAR = vcAR +colori*(-(g0L*muugl*sqrtS*(g1L*g2L*muu4*sqrtS*xc1 + 
     -       g1R*g2R*(Ml0*xc0 + (Ml0 + Ml2)*xc2))) + 
     -  g0R*(g1L*g2L*muu4*sqrtS*
     -      (Ml2*(xc1 + xc2) + Ml0*(xc0 + xc1 + xc2)) + 
     -     g1R*g2R*(xb0 + Ml0*(Ml0 + Ml2)*xc0 + 
     -        sqrtS**2*(muu4**2*xc1 + muu1**2*xc2))))


	  vcBL = vcBL + colori*g0R*(g1R*g2R*muu4*sqrtS*xc1 + g1L*g2L*(Ml2*xc2 + Ml0*(xc0 + xc2)))
	  vcBR = vcBR + colori*(-(g0L*(g1L*g2L*muu4*sqrtS*xc1 + g1R*g2R*(Ml2*xc2 + Ml0*(xc0 + xc2)))))
	
	  vcCL = vcCL + colori*(-(g1R*g2R*(-(g0R*muugl*sqrtS*xc1) + g0L*Ml2*(xc1 + xc2) + 
     -      g0L*Ml0*(xc0 + xc1 + xc2))))
     
	  vcCR = vcCL + colori*g1L*g2L*(-(g0L*muugl*sqrtS*xc1) + g0R*Ml2*(xc1 + xc2) + 
     -    g0R*Ml0*(xc0 + xc1 + xc2))
	  
	  vcDL = vcDL + colori*g0R*g1R*g2R*xc1
	  vcDR = vcDR + colori*(-(g0L*g1L*g2L*xc1)) 
	  
		  enddo
	     enddo
	    enddo
       enddo
	  enddo


! c-------------------------------Vertex3 ( Sfermion (j,l,n),Gluino, Gluon (M0,M1,M2))	
	  do j=1,2
	   do l=3,4
		do n=1,3

	
	  g0 = GS
	  g1L = GluLt1(j,l,ftt2,n,fgen2)
	  g1R = GluRt1(j,l,ftt2,n,fgen2)
	  g2 = VSfSf(1,isf1,j,itt1,l,igen1,n)
		
	  Ml0 = MSf(j,l,n)
	  Ml1 = MGl
	  Ml2 = 0d0
	
	  call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0
	
 	  vcAL = vcAL + colori*g0*g2*(g1R*muu4*sqrtS*
     -     (muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2)) - 
     -    g1L*(xb0 + Ml0**2*xc0 - 
     -       sqrtS*(muu4**2*sqrtS*(xc0 - xc2) + muu1**2*sqrtS*xc2 + 
     -          Ml1*muugl*(-xc0 + xc2))))
 	  
	  vcAR = vcAR + colori*g0*g2*(g1L*muu4*sqrtS*
     -     (muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2)) - 
     -    g1R*(xb0 + Ml0**2*xc0 - 
     -       sqrtS*(muu4**2*sqrtS*(xc0 - xc2) + muu1**2*sqrtS*xc2 + 
     -          Ml1*muugl*(-xc0 + xc2))))
	  
	  vcBL = vcBL + colori*g0*g2*(g1R*muu4*sqrtS*(xc0 + xc1 - xc2) + g1L*Ml1*(-xc0 + xc2))
     
	  vcBR = vcBR + colori*g0*g2*(g1L*muu4*sqrtS*(xc0 + xc1 - xc2) + g1R*Ml1*(-xc0 + xc2))
	
	
	  vcCL = vcCL + colori*g0*g1R*g2*(muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2))
	  vcCR = vcCL + colori*g0*g1L*g2*(muugl*sqrtS*(xc0 + xc1 - xc2) + Ml1*(-xc0 + xc1 + xc2))
	
	
	  vcDL = vcDL + colori*g0*g1R*g2*(xc0 + xc1 - xc2)
	  
	  vcDR = vcDR + colori*g0*g1L*g2*(xc0 + xc1 - xc2)

	    enddo
	   enddo
	  enddo


c-------------------------------Vertex4 ( Fermion (l,n), Gluon, Gluino (M0,M1,M2))	
! 
	  do l=3,4
	   do n=1,3
	
	  g0 = GS
	  g1 = BR(1,ftt2,l,fgen2,n)
	  g2L = GluLt1(isf1,itt1,l,igen1,n)
      g2R = GluRt1(isf1,itt1,l,igen1,n)
		
	  Ml0 = Mf(l,n)
	  Ml1 = 0d0
	  Ml2 = MGl
	
	  call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = -Nc/2d0

     
 	  vcAL = vcAL + colori*2*g0*g1*(g2R*Ml2*muu4*sqrtS*(xc1 + xc2) - 
     -    g2R*Ml0*(2*Ml2*xc0 + muugl*sqrtS*(xc0 + xc2)) + 
     -    g2L*(2*xb0 + 2*Ml0**2*xc0 - 
     -       Ml0*muu4*sqrtS*(xc0 + xc1 + xc2) + 
     -       sqrtS*(sqrtS*(muu2**2 + muu4**2 - tred)*xc1 + 
     -          Ml2*muugl*xc2 + 2*muu1**2*sqrtS*xc2)))
	   	  
	  vcAR = vcAR + colori*2*g0*g1*(g2L*Ml2*muu4*sqrtS*(xc1 + xc2) - 
     -    g2L*Ml0*(2*Ml2*xc0 + muugl*sqrtS*(xc0 + xc2)) + 
     -    g2R*(2*xb0 + 2*Ml0**2*xc0 - 
     -       Ml0*muu4*sqrtS*(xc0 + xc1 + xc2) + 
     -       sqrtS*(sqrtS*(muu2**2 + muu4**2 - tred)*xc1 + 
     -          Ml2*muugl*xc2 + 2*muu1**2*sqrtS*xc2)))
!      
	  vcBL = vcBL + colori*(-2*g0*g1*(-(g2L*Ml2*xc2) + g2R*Ml0*(xc0 + xc2)))
	  vcBR = vcBR + colori*(-2*g0*g1*(-(g2R*Ml2*xc2) + g2L*Ml0*(xc0 + xc2)))
	
	  vcCL = vcCL + colori*(-2*g0*g1*(-(g2R*Ml2*(xc1 + xc2)) + g2L*Ml0*(xc0 + xc1 + xc2)))
	  vcCR = vcCL + colori*(-2*g0*g1*(-(g2L*Ml2*(xc1 + xc2)) + g2R*Ml0*(xc0 + xc1 + xc2)))
	
	  vcDL = vcDL + colori* 0d0      
	
	  vcDR = vcDR + colori* 0d0
	
	   enddo
	  enddo

  	  vcAL = vcAL/(4d0*Pi)**2
  	  vcAR = vcAR/(4d0*Pi)**2
  	  vcBL = vcBL/(4d0*Pi)**2
  	  vcBR = vcBR/(4d0*Pi)**2
  	  vcCL = vcCL/(4d0*Pi)**2
  	  vcCR = vcCR/(4d0*Pi)**2
  	  vcDL = vcDL/(4d0*Pi)**2
  	  vcDR = vcDR/(4d0*Pi)**2
	  
  	
! ! ! 
! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! c -----------------------Amplitude square VC Gluino u-channel with gluino u-channel -----------
! ! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  colorf = Nc/2d0*CF
    
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1L*(g1Rc*g2Lc*muu3*
     -          (-(muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) + 
     -            muugl**2*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -             (-2*muu4*vcAR + 
     -               sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcBL) + 
     -            4*muu4*muugl*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g1Lc*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu4*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR)+
     -           g1Lc*g2Rc*(muu4**5*sqrtS*vcBL + 
     -            muu4**3*(-(muugl*vcAL) + 
     -               sqrtS*(-2 + 2*muu1**2 + muu2**2 + 3*muu3**2 + 
     -                  muugl**2 - 2*tred)*vcBL) + 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (-(muugl*vcAL) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muugl**2 - tred)*vcBL) - 
     -            muu4**4*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (muu1**2*(muu2**2 - tred) + 
     -               (muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred))
     -              *(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            muu4**2*(-1 + muu1**2 + 3*muu3**2 - 2*tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (-2*muu4*muugl*vcAL + 2*muu4**3*sqrtS*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu2**2 + muu3**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR))) + 
     -      g1R*(g1Lc*g2Lc*muu3*muugl*
     -          (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) - 2*muu4*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g1Rc*g2Lc*(-(muu1**2*muu2**2*vcAL) + muu3**2*vcAL - 
     -            muu2**2*muu3**2*vcAL - muu3**4*vcAL + 
     -            muu4**2*vcAL - muu1**2*muu4**2*vcAL - 
     -            3*muu3**2*muu4**2*vcAL - muu4**4*vcAL - tred*vcAL + 
     -            muu1**2*tred*vcAL + muu2**2*tred*vcAL + 
     -            2*muu3**2*tred*vcAL + 2*muu4**2*tred*vcAL - 
     -            tred**2*vcAL + muu4*muugl*vcAR - 
     -            muu1**2*muu4*muugl*vcAR - 
     -            2*muu3**2*muu4*muugl*vcAR - muu4**3*muugl*vcAR + 
     -            muu4*muugl*tred*vcAR + 
     -            2*muu1**2*muu2**2*muugl*sqrtS*vcBL - 
     -            2*muu3**2*muugl*sqrtS*vcBL + 
     -            2*muu2**2*muu3**2*muugl*sqrtS*vcBL + 
     -            2*muu3**4*muugl*sqrtS*vcBL - 
     -            2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muu1**2*muu4**2*muugl*sqrtS*vcBL + 
     -            6*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muu4**4*muugl*sqrtS*vcBL + 
     -            2*muugl*sqrtS*tred*vcBL - 
     -            2*muu1**2*muugl*sqrtS*tred*vcBL - 
     -            2*muu2**2*muugl*sqrtS*tred*vcBL - 
     -            4*muu3**2*muugl*sqrtS*tred*vcBL - 
     -            4*muu4**2*muugl*sqrtS*tred*vcBL + 
     -            2*muugl*sqrtS*tred**2*vcBL + 
     -            muu4*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g1Lc*g2Rc*muu3*
     -          (-2*muu4**3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            2*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            2*muu4**4*sqrtS*vcBR + 
     -            muu4**2*(-2*muugl*vcAR + 
     -               sqrtS*(-3 + 2*muu1**2 + 3*muu2**2 + 3*muu3**2 + 
     -                  2*muugl**2 - 3*tred)*vcBR) + 
     -            (-1 + muu2**2 + muu3**2 - tred)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muugl**2 - tred)*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(muu4**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu4**4*sqrtS*vcBR + 
     -            (-1 + muu3**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu4**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu3**2 + 
     -                  muugl**2 - tred)*vcBR)))))/
     -  (DenGlU*dconjg(DenGlU)))
! ! !       test = ststQQ_M2VertexT
! ! 	
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! c -----------------------Amplitude square VC Gluino u-channel with gluino t-channel -----------
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
	
	  colorf = -1d0/2d0*CF
	
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1R*(g1Rc*
     -          (-(g2Lc*muu3**2*muu4**2*vcAL) + g2Lc*tred*vcAL - 
     -            g2Lc*muu2**2*tred*vcAL - g2Lc*muu3**2*tred*vcAL - 
     -            g2Lc*muu4**2*tred*vcAL + g2Lc*tred**2*vcAL + 
     -            g2Lc*muu1**2*muu4*muugl*vcAR - 
     -            g2Lc*muu3**2*muu4*muugl*vcAR - 
     -            g2Rc*(-1 + muu3**2)*muugl**2*vcAR - 
     -            g2Rc*muu4**2*muugl**2*vcAR - 
     -            g2Lc*muu4*muugl*tred*vcAR + 
     -            2*g2Lc*muu3**2*muu4**2*muugl*sqrtS*vcBL - 
     -            2*g2Lc*muugl*sqrtS*tred*vcBL + 
     -            2*g2Lc*muu2**2*muugl*sqrtS*tred*vcBL + 
     -            2*g2Lc*muu3**2*muugl*sqrtS*tred*vcBL + 
     -            2*g2Lc*muu4**2*muugl*sqrtS*tred*vcBL - 
     -            2*g2Lc*muugl*sqrtS*tred**2*vcBL + 
     -            g2Lc*muu1**2*(muu2**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + 
     -            g2Rc*muu4**3*muugl*(-vcAL + 2*muugl*sqrtS*vcBL) + 
     -            g2Rc*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (-vcAL + 2*muugl*sqrtS*vcBL) + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*
     -             (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Lc*muu4*(-muu1**2 + muu3**2 + tred))*vcBR) + 
     -         g1Lc*muu3*(g2Lc*muugl*
     -             (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -                  vcAL) - 2*muu4*muugl*vcAR + 
     -               2*muugl*sqrtS*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL
     -                 + 2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            g2Rc*(muu1**2*muu4*vcAL + muu2**2*muu4*vcAL - 
     -               muu3**2*muu4*vcAL - muu4**3*vcAL + 
     -               muu2**2*muugl*vcAR - muu4**2*muugl*vcAR - 
     -               muugl*tred*vcAR - 
     -               2*muu1**2*muu4*muugl*sqrtS*vcBL - 
     -               2*muu2**2*muu4*muugl*sqrtS*vcBL + 
     -               2*muu3**2*muu4*muugl*sqrtS*vcBL + 
     -               2*muu4**3*muugl*sqrtS*vcBL + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 + muugl**2 - tred)*
     -                (-muu2**2 + muu4**2 + tred)*vcBR))) + 
     -      g1L*(g1Lc*g2Lc*muugl*
     -          (-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu4*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR)+
     -           g1Rc*g2Rc*muu3*muugl*
     -          (-2*muu4*muugl*vcAL + 2*muu4**3*sqrtS*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu2**2 + muu3**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1Rc*g2Lc*muu3*
     -          (-(muugl*tred*vcAL) - muu2**4*sqrtS*vcBL + 
     -            muu4**4*sqrtS*vcBL + 
     -            sqrtS*(-1 + muu1**2 + muu3**2 + muugl**2 - tred)*
     -             tred*vcBL + 
     -            muu4**2*(-(muugl*vcAL) + 
     -               (-1 + muu1**2 + muu3**2 + muugl**2)*sqrtS*vcBL) +
     -              (muu1 - muu3)*(muu1 + muu3)*muu4*
     -             (vcAR - 2*muugl*sqrtS*vcBR) - 
     -            muu4**3*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -            muu2**2*(muu4*vcAR - muugl**2*sqrtS*vcBL - 
     -               sqrtS*(-1 + muu1**2 + muu3**2 - 2*tred)*vcBL + 
     -               muugl*(vcAL - 2*muu4*sqrtS*vcBR))) + 
     -         g1Lc*g2Rc*(-(muu1**4*muu4*sqrtS*vcBL) + 
     -            muu3**4*muu4*sqrtS*vcBL + 
     -            tred*(-(muu4*muugl*vcAL) + muu4**3*sqrtS*vcBL + 
     -               muu4*sqrtS*(-1 + muu2**2 + muugl**2 - tred)*
     -                vcBL - muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -               (-1 + muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)
     -               ) + muu1**2*
     -             (muu4*muugl*vcAL - muu4**3*sqrtS*vcBL - 
     -               muu4*sqrtS*(-1 + muu2**2 + muugl**2 - 2*tred)*
     -                vcBL + 
     -               (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) + 
     -            muu3**2*(-(muu4*muugl*vcAL) + muu4**3*sqrtS*vcBL + 
     -               muu4*(-1 + muu2**2 + muugl**2)*sqrtS*vcBL - 
     -               muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -               tred*(vcAR - 2*muugl*sqrtS*vcBR))))))/
     -  (DenGlU*dconjg(DenGlT)))
! ! 
        do n=1,6
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! c -----------------------Amplitude square VC Gluino u-channel with Neutralino T-channel -----------
! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  colorf = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
! 

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1R*g3Lc*g4Rc*muu3*
     -       (muu1**2*muu4*vcAL + muu2**2*muu4*vcAL - 
     -         muu3**2*muu4*vcAL - muu4**3*vcAL + 
     -         muu2**2*muugl*vcAR - muu4**2*muugl*vcAR - 
     -         muugl*tred*vcAR - 2*muu1**2*muu4*muugl*sqrtS*vcBL - 
     -         2*muu2**2*muu4*muugl*sqrtS*vcBL + 
     -         2*muu3**2*muu4*muugl*sqrtS*vcBL + 
     -         2*muu4**3*muugl*sqrtS*vcBL + 
     -         sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -            muugl**2 - tred)*(-muu2**2 + muu4**2 + tred)*vcBR)+
     -        g1R*g3Rc*g4Lc*
     -       (-(muu1**4*muu4*sqrtS*vcBR) + muu3**4*muu4*sqrtS*vcBR + 
     -         muu3**2*(-(muu4**2*vcAL) - tred*vcAL - 
     -            muu4*muugl*vcAR + 2*muu4**2*muugl*sqrtS*vcBL + 
     -            2*muugl*sqrtS*tred*vcBL + 
     -            muu4*(-1 + muu2**2 + muu4**2 + muugl**2)*sqrtS*vcBR)
     -           + muu1**2*(muu4*muugl*vcAR + 
     -            (muu2**2 - tred)*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -            muu4*sqrtS*
     -             (-1 + muu2**2 + muu4**2 + muugl**2 - 2*tred)*vcBR)+
     -           tred*(-((-1 + muu2**2 + muu4**2 - tred)*vcAL) - 
     -            muu4*muugl*vcAR + 
     -            2*muugl*sqrtS*(-1 + muu2**2 + muu4**2 - tred)*
     -             vcBL + muu4*sqrtS*
     -             (-1 + muu2**2 + muu4**2 + muugl**2 - tred)*vcBR))+
     -        g1L*g3Rc*g4Lc*muu3*
     -       (-(muugl*tred*vcAL) - muu2**4*sqrtS*vcBL + 
     -         muu4**4*sqrtS*vcBL + 
     -         sqrtS*(-1 + muu1**2 + muu3**2 + muugl**2 - tred)*tred*
     -          vcBL + muu4**2*
     -          (-(muugl*vcAL) + 
     -            (-1 + muu1**2 + muu3**2 + muugl**2)*sqrtS*vcBL) + 
     -         (muu1 - muu3)*(muu1 + muu3)*muu4*
     -          (vcAR - 2*muugl*sqrtS*vcBR) - 
     -         muu4**3*(vcAR - 2*muugl*sqrtS*vcBR) + 
     -         muu2**2*(muu4*vcAR - muugl**2*sqrtS*vcBL - 
     -            sqrtS*(-1 + muu1**2 + muu3**2 - 2*tred)*vcBL + 
     -            muugl*(vcAL - 2*muu4*sqrtS*vcBR))) + 
     -      g1L*g3Lc*g4Rc*(-(muu1**4*muu4*sqrtS*vcBL) + 
     -         muu3**4*muu4*sqrtS*vcBL + 
     -         tred*(-(muu4*muugl*vcAL) + muu4**3*sqrtS*vcBL + 
     -            muu4*sqrtS*(-1 + muu2**2 + muugl**2 - tred)*vcBL - 
     -            muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)) +
     -           muu1**2*(muu4*muugl*vcAL - muu4**3*sqrtS*vcBL - 
     -            muu4*sqrtS*(-1 + muu2**2 + muugl**2 - 2*tred)*
     -             vcBL + (muu2**2 - tred)*(vcAR - 2*muugl*sqrtS*vcBR)
     -            ) + muu3**2*
     -          (-(muu4*muugl*vcAL) + muu4**3*sqrtS*vcBL + 
     -            muu4*(-1 + muu2**2 + muugl**2)*sqrtS*vcBL - 
     -            muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            tred*(vcAR - 2*muugl*sqrtS*vcBR))) + 
     -      (g1L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu4*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR)+
     -           g1R*g3Lc*g4Lc*muu3*
     -          (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) - 2*muu4*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g1L*g3Rc*g4Rc*muu3*
     -          (-2*muu4*muugl*vcAL + 2*muu4**3*sqrtS*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu2**2 + muu3**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1R*g3Rc*g4Rc*
     -          (-(muu4**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu4**4*sqrtS*vcBR + 
     -            (-1 + muu3**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu4**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu3**2 + 
     -                  muugl**2 - tred)*vcBR)))*muupt(n)))/
     -  (DenGlU*dconjg(DenGauT(n))))

! ! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! c -----------------------Amplitude square VC Gluino u-channel with Neutralino U-channel -----------
! ! c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!       
! ! 	  
	  colorf= 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1L*g3Rc*g4Lc*muu3*
     -       (-(muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -            vcAL) + muugl**2*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -         (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-2*muu4*vcAR + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcBL) + 4*muu4*muugl*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g1R*g3Rc*g4Lc*
     -       (-(muu1**2*muu2**2*vcAL) + muu3**2*vcAL - 
     -         muu2**2*muu3**2*vcAL - muu3**4*vcAL + muu4**2*vcAL - 
     -         muu1**2*muu4**2*vcAL - 3*muu3**2*muu4**2*vcAL - 
     -         muu4**4*vcAL - tred*vcAL + muu1**2*tred*vcAL + 
     -         muu2**2*tred*vcAL + 2*muu3**2*tred*vcAL + 
     -         2*muu4**2*tred*vcAL - tred**2*vcAL + muu4*muugl*vcAR - 
     -         muu1**2*muu4*muugl*vcAR - 2*muu3**2*muu4*muugl*vcAR - 
     -         muu4**3*muugl*vcAR + muu4*muugl*tred*vcAR + 
     -         2*muu1**2*muu2**2*muugl*sqrtS*vcBL - 
     -         2*muu3**2*muugl*sqrtS*vcBL + 
     -         2*muu2**2*muu3**2*muugl*sqrtS*vcBL + 
     -         2*muu3**4*muugl*sqrtS*vcBL - 
     -         2*muu4**2*muugl*sqrtS*vcBL + 
     -         2*muu1**2*muu4**2*muugl*sqrtS*vcBL + 
     -         6*muu3**2*muu4**2*muugl*sqrtS*vcBL + 
     -         2*muu4**4*muugl*sqrtS*vcBL + 2*muugl*sqrtS*tred*vcBL - 
     -         2*muu1**2*muugl*sqrtS*tred*vcBL - 
     -         2*muu2**2*muugl*sqrtS*tred*vcBL - 
     -         4*muu3**2*muugl*sqrtS*tred*vcBL - 
     -         4*muu4**2*muugl*sqrtS*tred*vcBL + 
     -         2*muugl*sqrtS*tred**2*vcBL + 
     -         muu4*sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -            muugl**2 - tred)*vcBR) + 
     -      g1L*g3Lc*g4Rc*(muu4**5*sqrtS*vcBL + 
     -         muu4**3*(-(muugl*vcAL) + 
     -            sqrtS*(-2 + 2*muu1**2 + muu2**2 + 3*muu3**2 + 
     -               muugl**2 - 2*tred)*vcBL) + 
     -         muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -          (-(muugl*vcAL) + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -               muugl**2 - tred)*vcBL) - 
     -         muu4**4*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -         (muu1**2*(muu2**2 - tred) + 
     -            (muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred))*
     -          (vcAR - 2*muugl*sqrtS*vcBR) - 
     -         muu4**2*(-1 + muu1**2 + 3*muu3**2 - 2*tred)*
     -          (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -      g1R*g3Lc*g4Rc*muu3*
     -       (-2*muu4**3*(vcAL - 2*muugl*sqrtS*vcBL) - 
     -         2*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 - tred)*
     -          (vcAL - 2*muugl*sqrtS*vcBL) + 2*muu4**4*sqrtS*vcBR + 
     -         muu4**2*(-2*muugl*vcAR + 
     -            sqrtS*(-3 + 2*muu1**2 + 3*muu2**2 + 3*muu3**2 + 
     -               2*muugl**2 - 3*tred)*vcBR) + 
     -         (-1 + muu2**2 + muu3**2 - tred)*
     -          (-(muugl*vcAR) + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -               muugl**2 - tred)*vcBR)) + 
     -      (g1L*g3Lc*g4Lc*(-((-1 + muu3**2 + muu4**2)*muugl*vcAL) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcAR + (-1 + muu3**2 + muu4**2)*muugl**2*sqrtS*
     -             vcBL + (-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            2*muu4*muugl*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR)+ 
     -          g1R*g3Lc*g4Lc*muu3*
     -          (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) - 2*muu4*muugl*vcAR + 
     -            2*muugl*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 + 
     -               muugl**2 - tred)*vcBR) + 
     -         g1L*g3Rc*g4Rc*muu3*
     -          (-2*muu4*muugl*vcAL + 2*muu4**3*sqrtS*vcBL + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muugl**2 - 
     -               tred)*vcBL - 
     -            2*muu4**2*(vcAR - 2*muugl*sqrtS*vcBR) - 
     -            (-1 + muu2**2 + muu3**2 - tred)*
     -             (vcAR - 2*muugl*sqrtS*vcBR)) + 
     -         g1R*g3Rc*g4Rc*
     -          (-(muu4**3*(vcAL - 2*muugl*sqrtS*vcBL)) - 
     -            muu4*(-1 + muu1**2 + 2*muu3**2 - tred)*
     -             (vcAL - 2*muugl*sqrtS*vcBL) + muu4**4*sqrtS*vcBR + 
     -            (-1 + muu3**2)*
     -             (-(muugl*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muugl**2 - tred)*vcBR) + 
     -            muu4**2*(-(muugl*vcAR) + 
     -               sqrtS*(-2 + muu1**2 + muu2**2 + 2*muu3**2 + 
     -                  muugl**2 - tred)*vcBR)))*muupu(n)))/
     -  (DenGlU*dconjg(DenGauU(n))))     
! 
! ! 
! 
        enddo

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

c****************************************************
c************ Loop external squark (ifs1,itt1,igen1) & quark (ftt1,fgen1) with Neutralino in t-channel upper vertex ***********************
c****************************************************	
	
	  do neu=1,6
	
	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0
	
c-------------------------------Vertex1 ( Gluon , fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
	  do i=1,2
	   do k=3,4
	    do l=3,4
        m = 3
        n = 3
 
	  g0L = AL(neu,i,k,l,m,n)
	  g0R = AR(neu,i,k,l,m,n)
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2 = VSfSf(1,isf1,i,itt1,k,igen1,m)

	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)	
	
	  colori = CF

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-tman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu3*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))
     
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-tman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu3*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))

     
	  vcBL = vcBL + colori*g1*g2*g0R*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	 
	  vcCL = vcCL + colori*g1*g2*(g0R*MGau(neu)*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGau(neu)*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)

	    enddo
	   enddo
	  enddo

c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
       do k=3,4
	    do l=3,4
	    m = 3
	    n = 3
!  
	  g0L = dconjg(AR(neu,j,l,k,n,m))
	  g0R = dconjg(AL(neu,j,l,k,n,m))
	  g1L = GluR(j,l,ftt1,n,fgen1)
	  g1R = GluL(j,l,ftt1,n,fgen1)
	  g2L = GluR(isf1,itt1,k,igen1,m)
	  g2R = GluL(isf1,itt1,k,igen1,m)
	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu3**2*sqrtS**2,tman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	  
	
	  colori = CF
! 	
 	  vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGau(neu)*muu3*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGau(neu)*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGau(neu)*(xc2+xc0)*g1L*g2L*g0R))
     
	  vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGau(neu)*muu3*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGau(neu)*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGau(neu)*(xc2+xc0)*g1R*g2R*g0L))


	  vcBL = vcBL - colori* (muu3*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	  vcBR = vcBR - colori* (muu3*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	  vcCL = vcCL - colori* (MGau(neu)*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	  vcCR = vcCL - colori* (MGau(neu)*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	  vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	  vcDR = vcDR + colori* xc1*g1L*g2R*g0L


	     enddo
	    enddo
	  enddo	



! 	
	vcAL = vcAL/(4d0*Pi)**2
 	vcAR = vcAR/(4d0*Pi)**2
 	
 	vcBL = vcBL/(4d0*Pi)**2
 	vcBR = vcBR/(4d0*Pi)**2
 	
 	vcCL = vcCL/(4d0*Pi)**2
 	vcCR = vcCR/(4d0*Pi)**2
 	
	vcDL = vcDL/(4d0*Pi)**2
 	vcDR = vcDR/(4d0*Pi)**2
	
c****************************************************
c************ Amplitudenquadrate*********************
c****************************************************

c-----------------------T-channel Neutralino squared------------------------------------------------------- 
	  do neuc=1,6
	  
      n=neuc
      k=neu
      
      colorf = Nc**2
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3Lc*g4Rc*
     -       (-(g4R*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAL) + 
     -         g4L*muu4*tred*
     -          (-2*muu3*vcAR + 
     -            sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) + 
     -         g4R*muu3*sqrtS*tred*(-muu2**2 + muu4**2 + tred)*vcBR) +
     -       g3Rc*g4Lc*(-(g4L*
     -            ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAR) + 
     -         g4L*muu3*sqrtS*tred*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -         g4R*muu4*tred*
     -          (-2*muu3*vcAL + 
     -            sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)) + 
     -      muuqt(neu)*(g3Rc*g4Lc*
     -          (g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -            g4R*muu4*(muu1**2*vcAR - (muu3**2 + tred)*vcAR + 
     -               4*muu3*sqrtS*tred*vcBL) + 
     -            2*g4L*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBR) + 
     -         g3Lc*g4Rc*(g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -            2*g4R*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBL + 
     -            g4L*muu4*(muu1**2*vcAL - (muu3**2 + tred)*vcAL + 
     -               4*muu3*sqrtS*tred*vcBR)) + 
     -         sqrtS*(g3Rc*g4Lc*
     -             (g4L*muu3*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -               g4R*muu4*(-muu1**2 + muu3**2 + tred)*vcBR) + 
     -            g3Lc*g4Rc*(g4L*muu4*(-muu1**2 + muu3**2 + tred)*
     -                vcBL + g4R*muu3*(-muu2**2 + muu4**2 + tred)*vcBR
     -               ))*muuqt(neu)) + 
     -      muupt(neuc)*(g3Rc*g4Rc*
     -          (g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -            g4L*muu4*(muu1**2*vcAR - (muu3**2 + tred)*vcAR + 
     -               2*muu3*sqrtS*tred*vcBL) + 
     -            g4R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -         g3Lc*g4Lc*(g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -            g4L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            g4R*muu4*(muu1**2*vcAL - (muu3**2 + tred)*vcAL + 
     -               2*muu3*sqrtS*tred*vcBR)) + 
     -         muuqt(neu)*(g3Lc*g4Lc*
     -             (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -               2*g4R*muu4*
     -                (-(muu3*vcAR) + 
     -                  sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) + 
     -               2*g4L*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR
     -               ) + g3Rc*g4Rc*
     -             (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -               2*g4R*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*
     -                vcBL + 
     -               2*g4L*muu4*
     -                (-(muu3*vcAL) + 
     -                  sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)) + 
     -            sqrtS*(g3Lc*g4Lc*
     -                (g4L*(-1 + muu3**2 + muu4**2)*vcBL + 
     -                  2*g4R*muu3*muu4*vcBR) + 
     -               g3Rc*g4Rc*
     -                (2*g4L*muu3*muu4*vcBL + 
     -                  g4R*(-1 + muu3**2 + muu4**2)*vcBR))*muuqt(neu)))
     -      ))/(dconjg(DenGauT(neuc))*DenGauT(neu)))
!       
! ! c-----------------------T-channel correction Neutralino with u-channel neutralino------------------------------------------------------- 
! 
      colorf = Nc
	
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
! 	
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3Lc*g4Rc*
     -       (g4R*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAL + 
     -         g4L*muu4*(muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -             vcAR + sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*tred*
     -             vcBL) + g4R*muu3*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*vcBR) +
     -         g3Rc*g4Lc*(g4L*
     -          (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAR + 
     -         g4L*muu3*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*vcBL +
     -           g4R*muu4*(muu3*
     -             (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             tred*vcBR)) + 
     -      muuqt(neu)*(g3Rc*g4Lc*
     -          (-(g4L*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                 tred)*vcAL) + 
     -            g4R*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu3*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -                sqrtS*vcBL) + 
     -            2*g4L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g3Lc*g4Rc*(-(g4R*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR)
     -              + 2*g4R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g4L*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAL) + 
     -               2*muu3*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -                sqrtS*vcBR)) + 
     -         sqrtS*(g3Rc*g4Lc*
     -             (g4L*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL
     -                 + g4R*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR
     -               ) + g3Lc*g4Rc*
     -             (g4L*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL
     -                 + g4R*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR
     -               ))*muuqt(neu)) + 
     -      muupu(neuc)*(g3Rc*g4Rc*
     -          (g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -            g4L*muu4*(muu1**2*vcAR - (muu3**2 + tred)*vcAR + 
     -               2*muu3*sqrtS*tred*vcBL) + 
     -            g4R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -         g3Lc*g4Lc*(g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -            g4L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            g4R*muu4*(muu1**2*vcAL - (muu3**2 + tred)*vcAL + 
     -               2*muu3*sqrtS*tred*vcBR)) + 
     -         muuqt(neu)*(g3Lc*g4Lc*
     -             (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -               2*g4R*muu4*
     -                (-(muu3*vcAR) + 
     -                  sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) + 
     -               2*g4L*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR
     -               ) + g3Rc*g4Rc*
     -             (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -               2*g4R*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*
     -                vcBL + 
     -               2*g4L*muu4*
     -                (-(muu3*vcAL) + 
     -                  sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)) + 
     -            sqrtS*(g3Lc*g4Lc*
     -                (g4L*(-1 + muu3**2 + muu4**2)*vcBL + 
     -                  2*g4R*muu3*muu4*vcBR) + 
     -               g3Rc*g4Rc*
     -                (2*g4L*muu3*muu4*vcBL + 
     -                  g4R*(-1 + muu3**2 + muu4**2)*vcBR))*muuqt(neu)))
     -      ))/(dconjg(DenGauU(neuc))*DenGauT(neu)))
      
      
      enddo !Index neuc
      
c-----------------------T-channel correction Neutralino with t-channel gluino------------------------------------------------------- 

      colorf= 0d0
	  
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -         g4L*muu4*(muu1**2*vcAR - (muu3**2 + tred)*vcAR + 
     -            2*muu3*sqrtS*tred*vcBL) + 
     -         g4R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -      g1Lc*g2Rc*(-(g4R*
     -            ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAL) + 
     -         g4L*muu4*tred*
     -          (-2*muu3*vcAR + 
     -            sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) + 
     -         g4R*muu3*sqrtS*tred*(-muu2**2 + muu4**2 + tred)*vcBR) +
     -       g1Lc*g2Lc*muugl*
     -       (g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -         g4L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -         g4R*muu4*(muu1**2*vcAL - (muu3**2 + tred)*vcAL + 
     -            2*muu3*sqrtS*tred*vcBR)) + 
     -      g1Rc*g2Lc*(-(g4L*
     -            ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAR) + 
     -         g4L*muu3*sqrtS*tred*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -         g4R*muu4*tred*
     -          (-2*muu3*vcAL + 
     -            sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)) - 
     -      (g1Lc*g2Lc*muugl*
     -          (g4L*(-1 + muu3**2 + muu4**2)*vcAL + 
     -            2*g4R*muu4*
     -             (muu3*vcAR - 
     -               sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) - 
     -            2*g4L*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR) +
     -           g1Rc*g2Lc*(g4L*muu3*(-muu2**2 + muu4**2 + tred)*
     -             vcAL + g4R*muu4*
     -             (-(muu1**2*vcAR) + (muu3**2 + tred)*vcAR - 
     -               4*muu3*sqrtS*tred*vcBL) - 
     -            2*g4L*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBR) + 
     -         g1Lc*g2Rc*(g4R*muu3*(-muu2**2 + muu4**2 + tred)*vcAR - 
     -            2*g4R*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBL + 
     -            g4L*muu4*(-(muu1**2*vcAL) + (muu3**2 + tred)*vcAL - 
     -               4*muu3*sqrtS*tred*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (g4R*(-1 + muu3**2 + muu4**2)*vcAR - 
     -            2*g4R*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -            2*g4L*muu4*
     -             (muu3*vcAL - 
     -               sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)))*
     -       muuqt(neu) + sqrtS*
     -       (g4L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -            g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Lc*g2Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g1Rc*g2Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL + 
     -         g4R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -            g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Rc*g2Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g1Lc*g2Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)*
     -       muuqt(neu)**2))/(dconjg(DenGlT)*DenGauT(neu)))
      
      
      
      
c-----------------------T-channel correction Neutralino with u-channel gluino------------------------------------------------------- 

      colorf = CF*Nc
	  
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  

	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -         g4L*muu4*(muu1**2*vcAR - (muu3**2 + tred)*vcAR + 
     -            2*muu3*sqrtS*tred*vcBL) + 
     -         g4R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -      g1Lc*g2Rc*(g4R*(muu1**2*muu2**2 - muu3**2*muu4**2 + 
     -            tred - (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             tred + tred**2)*vcAL + 
     -         g4L*muu4*(muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -             vcAR + sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*tred*
     -             vcBL) + g4R*muu3*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*vcBR)
     -        + g1Lc*g2Lc*muugl*
     -       (g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -         g4L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -         g4R*muu4*(muu1**2*vcAL - (muu3**2 + tred)*vcAL + 
     -            2*muu3*sqrtS*tred*vcBR)) + 
     -      g1Rc*g2Lc*(g4L*(muu1**2*muu2**2 - muu3**2*muu4**2 + 
     -            tred - (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             tred + tred**2)*vcAR + 
     -         g4L*muu3*sqrtS*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*vcBL+
     -           g4R*muu4*(muu3*
     -             (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             tred*vcBR)) - 
     -      (g1Lc*g2Lc*muugl*
     -          (g4L*(-1 + muu3**2 + muu4**2)*vcAL + 
     -            2*g4R*muu4*
     -             (muu3*vcAR - 
     -               sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL) - 
     -            2*g4L*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR) +
     -           g1Rc*g2Lc*(g4L*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAL + 
     -            g4R*muu4*((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcAR + 
     -               2*muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -                sqrtS*vcBL) - 
     -            2*g4L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g1Lc*g2Rc*(g4R*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR - 
     -            2*g4R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g4L*muu4*((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcAL + 
     -               2*muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -                sqrtS*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (g4R*(-1 + muu3**2 + muu4**2)*vcAR - 
     -            2*g4R*muu3*sqrtS*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -            2*g4L*muu4*
     -             (muu3*vcAL - 
     -               sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)))*
     -       muuqt(k) + sqrtS*
     -       (g4L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -            g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Lc*g2Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g1Rc*g2Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcBL +
     -           g4R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -            g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Rc*g2Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g1Lc*g2Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcBR)*
     -       muuqt(k)**2))/(dconjg(DenGlU)*DenGauT(k)))      
      
      
      enddo !Index neu

      
c****************************************************
c************ lower Loop external squark (ifs2,itt2,igen2) & quark (ftt2,fgen2) with Neutralino in t-channel ***********************
c****************************************************
       
	  do neu = 1,6
	  
	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0
	  	
c-------------------------------Vertex1 ( Gluon , Fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
      do i=1,2
	   do k=3,4
	    do l=3,4
        m = 3
        n = 3

	  g0L = AL(neu,i,k,l,m,n)
	  g0R = AR(neu,i,k,l,m,n)
	  g1 = BR(1,ftt2,l,fgen2,n)
	  g2 = VSfSf(1,isf2,i,itt2,k,igen2,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)	
	
	  colori = CF

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-tman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu4*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))
     
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-tman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu4*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))
!      
!      
	  vcBL = vcBL + colori*g1*g2*g0R*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	 
	  vcCL = vcCL + colori*g1*g2*(g0R*MGau(neu)*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGau(neu)*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)

	    enddo
	   enddo
	  enddo




c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	 do j=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
		  do n=1,3
		  
      g0L = dconjg(ARt1(neu,j,l,k,n,m))
	  g0R = dconjg(ALt1(neu,j,l,k,n,m))
	  g1L = GluR(j,l,ftt2,n,fgen2)
	  g1R = GluL(j,l,ftt2,n,fgen2)
	  g2L = GluR(isf2,itt2,k,igen2,m)
	  g2R = GluL(isf2,itt2,k,igen2,m)
	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu4**2*sqrtS**2,tman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(tman,Ml1**2,Ml2**2)
	
	  colori = CF
! 	
 	  vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu4**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGau(neu)*muu4*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu4*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGau(neu)*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu4*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGau(neu)*(xc2+xc0)*g1L*g2L*g0R))
     
	  vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu4**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGau(neu)*muu4*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu4*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGau(neu)*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu4*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGau(neu)*(xc2+xc0)*g1R*g2R*g0L))


	  vcBL = vcBL - colori* (muu4*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	  vcBR = vcBR - colori* (muu4*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	  vcCL = vcCL - colori* (MGau(neu)*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	  vcCR = vcCL - colori* (MGau(neu)*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	  vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	  vcDR = vcDR + colori* xc1*g1L*g2R*g0L
	  

		  enddo
	     enddo
	    enddo
       enddo
	  enddo
! 	  

! 		 

  	  vcAL = vcAL/(4d0*Pi)**2
  	  vcAR = vcAR/(4d0*Pi)**2
  	  vcBL = vcBL/(4d0*Pi)**2
  	  vcBR = vcBR/(4d0*Pi)**2
  	  vcCL = vcCL/(4d0*Pi)**2
  	  vcCR = vcCR/(4d0*Pi)**2
  	  vcDL = vcDL/(4d0*Pi)**2
  	  vcDR = vcDR/(4d0*Pi)**2

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c -----------------------Amplitude square VC Neutralino t-channel with Neutralino T-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do neuc = 1,6

      
	  colorf = Nc**2
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g3R*(-2*g3Lc*g4Rc*muu3*muu4*tred + 
     -         g3Rc*g4Lc*(-((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4)) + (-1 + muu1**2 + muu2**2)*tred - 
     -            tred**2))*vcAL + 
     -      g3L*(-2*g3Rc*g4Lc*muu3*muu4*tred + 
     -         g3Lc*g4Rc*(-((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4)) + (-1 + muu1**2 + muu2**2)*tred - 
     -            tred**2))*vcAR + 
     -      g3L*sqrtS*tred*(g3Lc*g4Rc*muu4*
     -          (-muu1**2 + muu3**2 + tred) + 
     -         g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL + 
     -      g3R*sqrtS*tred*(g3Rc*g4Lc*muu4*
     -          (-muu1**2 + muu3**2 + tred) + 
     -         g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR + 
     -      muuqt(neu)*(g3L*(g3Lc*g4Rc*muu4*
     -             (muu1**2 - muu3**2 - tred) + 
     -            g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*vcAL + 
     -         g3Lc*g3R*g4Rc*muu3*
     -          (muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -            4*muu4*sqrtS*tred*vcBL) + 
     -         g3R*g3Rc*g4Lc*
     -          (muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            2*sqrtS*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBL) + 
     -         2*g3L*sqrtS*(2*g3Rc*g4Lc*muu3*muu4*tred + 
     -            g3Lc*g4Rc*((muu1 - muu3)*(muu1 + muu3)*
     -                (muu2 - muu4)*(muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2))*vcBR + 
     -         sqrtS*(g3L*(g3Lc*g4Rc*muu4*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -              g3R*(g3Rc*g4Lc*muu4*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)
     -           *muuqt(neu)) + 
     -      muupt(neuc)*(g3Rc*g4Rc*
     -          (g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAL + 
     -            g3L*muu3*(muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -               2*muu4*sqrtS*tred*vcBL) + 
     -            g3R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -         g3Lc*g4Lc*(g3L*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            g3L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            g3R*muu3*(muu2**2*vcAL - (muu4**2 + tred)*vcAL + 
     -               2*muu4*sqrtS*tred*vcBR)) + 
     -         muuqt(neu)*(-(g3L*
     -               (2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcAL) - 
     -            g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcAR + 
     -            2*g3R*sqrtS*
     -             (g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -              2*g3L*sqrtS*
     -             (g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR +
     -              sqrtS*(g3L*
     -                (2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcBL + 
     -               g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcBR)*
     -             muuqt(neu)))))/(dconjg(DenGauT(neuc))*DenGauT(neu)))

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c -----------------------Amplitude square VC Nuetralino t-channel with neutralino U-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      colorf = Nc
	
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g3R*(g3Lc*g4Rc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Rc*g4Lc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2))*vcAL + 
     -      g3L*(g3Rc*g4Lc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Lc*g4Rc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2))*vcAR + 
     -      g3L*sqrtS*(g3Lc*g4Rc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Rc*g4Lc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*tred*vcBL
     -        + g3R*sqrtS*(g3Rc*g4Lc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Lc*g4Rc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*tred*vcBR
     -        + muuqt(neu)*(-(g3L*
     -            (g3Lc*g4Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3Rc*g4Lc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcAL
     -            ) - g3R*(g3Rc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Lc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcAR+
     -           2*g3R*sqrtS*
     -          (g3Lc*g4Rc*muu3*muu4*
     -             (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g3Rc*g4Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))*vcBL + 
     -         2*g3L*sqrtS*(g3Rc*g4Lc*muu3*muu4*
     -             (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g3Lc*g4Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))*vcBR + 
     -         sqrtS*(g3L*(g3Lc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Rc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + g3R*
     -             (g3Rc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Lc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR)*muuqt(neu)) + 
     -      muupu(neuc)*(g3Rc*g4Rc*
     -          (g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAL + 
     -            g3L*muu3*(muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -               2*muu4*sqrtS*tred*vcBL) + 
     -            g3R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -         g3Lc*g4Lc*(g3L*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            g3L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -            g3R*muu3*(muu2**2*vcAL - (muu4**2 + tred)*vcAL + 
     -               2*muu4*sqrtS*tred*vcBR)) + 
     -         muuqt(neu)*(-(g3L*
     -               (2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcAL) - 
     -            g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcAR + 
     -            2*g3R*sqrtS*
     -             (g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -              2*g3L*sqrtS*
     -             (g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR +
     -              sqrtS*(g3L*
     -                (2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcBL + 
     -               g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcBR)*
     -             muuqt(neu)))))/(dconjg(DenGauU(neuc))*DenGauT(neu)))
        enddo !end neuc
 
     
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c -----------------------Amplitude square VC Neutralino t-channel with Gluino T-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  colorf = 0d0
	  
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  
	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	
	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)


	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal((sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAL + 
     -         g3L*muu3*(muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -            2*muu4*sqrtS*tred*vcBL) + 
     -         g3R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -      g1Rc*g2Lc*(-(g3R*
     -            ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAL) + 
     -         g3L*muu3*tred*
     -          (-2*muu4*vcAR + 
     -            sqrtS*(-muu2**2 + muu4**2 + tred)*vcBL) + 
     -         g3R*muu4*sqrtS*tred*(-muu1**2 + muu3**2 + tred)*vcBR)+
     -        g1Lc*g2Lc*muugl*
     -       (g3L*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -         g3L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -         g3R*muu3*(muu2**2*vcAL - (muu4**2 + tred)*vcAL + 
     -            2*muu4*sqrtS*tred*vcBR)) + 
     -      g1Lc*g2Rc*(-(g3L*
     -            ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -               (muu2 + muu4) + tred - 
     -              (muu1**2 + muu2**2)*tred + tred**2)*vcAR) + 
     -         g3L*muu4*sqrtS*tred*(-muu1**2 + muu3**2 + tred)*vcBL + 
     -         g3R*muu3*tred*
     -          (-2*muu4*vcAL + 
     -            sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)) + 
     -      muuqt(neu)*(g1Lc*g2Lc*muugl*
     -          (-(g3L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -            2*g3R*muu3*
     -             (-(muu4*vcAR) + 
     -               sqrtS*(-muu2**2 + muu4**2 + tred)*vcBL) + 
     -            2*g3L*muu4*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR)+
     -           g1Lc*g2Rc*(g3L*muu4*(muu1**2 - muu3**2 - tred)*
     -             vcAL + g3R*muu3*
     -             (muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -               4*muu4*sqrtS*tred*vcBL) + 
     -            2*g3L*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBR) + 
     -         g1Rc*g2Lc*(g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            2*g3R*sqrtS*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)*vcBL + 
     -            g3L*muu3*(muu2**2*vcAL - (muu4**2 + tred)*vcAL + 
     -               4*muu4*sqrtS*tred*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(g3R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -            2*g3R*muu4*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL + 
     -            2*g3L*muu3*
     -             (-(muu4*vcAL) + 
     -               sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)) + 
     -         sqrtS*(g3L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -               g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Lc*g2Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Rc*g2Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL+
     -             g3R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -               g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Rc*g2Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Lc*g2Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)
     -           *muuqt(neu))))/(dconjg(DenGlT)*DenGauT(neu)))

     
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c -----------------------Amplitude square VC Neutralino t-channel with Gluino U-channel -----------
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      colorf = CF*Nc
	  
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))

	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

	  g3L = ALt1(neu,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(neu,isf1,itt1,ftt1,igen1,fgen1)

	  g4L = ALt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  g4R = ARt1(neu,isf2,itt2,ftt2,igen2,fgen2)
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*
     -		dreal( (sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAL + 
     -         g3L*muu3*(muu2**2*vcAR - (muu4**2 + tred)*vcAR + 
     -            2*muu4*sqrtS*tred*vcBL) + 
     -         g3R*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBR) + 
     -      g1Rc*g2Lc*(g3R*(muu1**2*muu2**2 - muu3**2*muu4**2 + 
     -            tred - (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             tred + tred**2)*vcAL + 
     -         g3L*muu3*(muu4*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -             vcAR + sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred*
     -             vcBL) + g3R*muu4*sqrtS*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*tred*vcBR)
     -        + g1Lc*g2Lc*muugl*
     -       (g3L*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -         g3L*(-1 + muu3**2 + muu4**2)*sqrtS*tred*vcBL + 
     -         g3R*muu3*(muu2**2*vcAL - (muu4**2 + tred)*vcAL + 
     -            2*muu4*sqrtS*tred*vcBR)) + 
     -      g1Lc*g2Rc*(g3L*(muu1**2*muu2**2 - muu3**2*muu4**2 + 
     -            tred - (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             tred + tred**2)*vcAR + 
     -         g3L*muu4*sqrtS*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*tred*vcBL +
     -           g3R*muu3*(muu4*
     -             (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             tred*vcBR)) - 
     -      (g1Lc*g2Lc*muugl*
     -          (g3L*(-1 + muu3**2 + muu4**2)*vcAL + 
     -            2*g3R*muu3*
     -             (muu4*vcAR - 
     -               sqrtS*(-muu2**2 + muu4**2 + tred)*vcBL) - 
     -            2*g3L*muu4*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBR) +
     -           g1Lc*g2Rc*(g3L*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAL + 
     -            g3R*muu3*((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcAR + 
     -               2*muu4*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -                sqrtS*vcBL) - 
     -            2*g3L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g1Rc*g2Lc*(g3R*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR - 
     -            2*g3R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g3L*muu3*((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcAL + 
     -               2*muu4*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -                sqrtS*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (g3R*(-1 + muu3**2 + muu4**2)*vcAR - 
     -            2*g3R*muu4*sqrtS*(-muu1**2 + muu3**2 + tred)*vcBL + 
     -            2*g3L*muu3*
     -             (muu4*vcAL - 
     -               sqrtS*(-muu2**2 + muu4**2 + tred)*vcBR)))*
     -       muuqt(neu) + sqrtS*
     -       (g3L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -            g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Lc*g2Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g1Rc*g2Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcBL +
     -           g3R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -            g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g1Rc*g2Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g1Lc*g2Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcBR)*
     -       muuqt(neu)**2))/(dconjg(DenGlU)*DenGauT(neu)))    

        enddo

        
c****************************************************
c************ Loop external squark (ifs2,itt2,igen2) & quark (ftt1,fgen1) with Neutralino in u-channel upper vertex ***********************
c****************************************************	

      do neu=1,6

	  vcAR = 0d0
	  vcAL = 0d0
	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0
	

c-------------------------------Vertex1 ( Gluon , fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
	  do i=1,2
	   do k=3,4
	    do l=3,4
	     do m=1,3
          do n=1,3

          
	  g0L = ALt1(neu,i,k,l,m,n)
	  g0R = ARt1(neu,i,k,l,m,n)
	  g1 = BR(1,ftt1,l,fgen1,n)
	  g2 = VSfSf(1,isf2,i,itt2,k,igen2,m)
	
	  Ml0 = 0d0
	  Ml1 = Mf(l,n)
	  Ml2 = MSf(i,k,m)
	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)	
	
	  colori = CF

	  vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-uman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu3*sqrtS*muuqu(neu)*sqrtS*(2*xc1+xc2+2*xc0)-muuqu(neu)*sqrtS*Ml1*(xc2+2*xc0)))
	     
	  vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu2**2*sqrtS**2*xc2 + 
     -     muu3**2*sqrtS**2*xc1+(-uman+muu2**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu3*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu3*sqrtS*muuqu(neu)*sqrtS*(2*xc1+xc2+2*xc0)-muuqu(neu)*sqrtS*Ml1*(xc2+2*xc0)))
!      
!      
	  vcBL = vcBL + colori*g1*g2*g0R*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcBR = vcBR + colori*g1*g2*g0L*(muu3*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	  vcCL = vcCL + colori*g1*g2*(g0R*MGau(neu)*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	  vcCR = vcCR + colori*g1*g2*(g0L*MGau(neu)*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	  vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	  vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)

		  enddo
	     enddo
	    enddo
	   enddo
	  enddo

! ! c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
         do m=1,3
		  do n=1,3


	  g0L = dconjg(ARt1(neu,j,l,k,n,m))
	  g0R = dconjg(ALt1(neu,j,l,k,n,m))
	  g1L = GluLt1(j,l,ftt1,n,fgen1)
	  g1R = GluRt1(j,l,ftt1,n,fgen1)
	  g2L = GluLt1(isf2,itt2,k,igen2,m)
	  g2R = GluRt1(isf2,itt2,k,igen2,m)
	
	  Ml0 = MGl
	  Ml1 = MSf(j,l,n)
	  Ml2 = Mf(k,m)
	
	  call myCget(muu3**2*sqrtS**2,uman,muu2**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	  call myBget(uman,Ml1**2,Ml2**2)
	
	  colori = CF
! 	
 	  vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGau(neu)*muu3*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGau(neu)*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGau(neu)*(xc2+xc0)*g1L*g2L*g0R))
    
     
	  vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu3**2*sqrtS**2*xc1 + muu2**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGau(neu)*muu3*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu3*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGau(neu)*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu3*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGau(neu)*(xc2+xc0)*g1R*g2R*g0L))


	  vcBL = vcBL - colori* (muu3*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	  vcBR = vcBR - colori* (muu3*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	  vcCL = vcCL - colori* (MGau(neu)*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	  vcCR = vcCL - colori* (MGau(neu)*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	  vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	  vcDR = vcDR + colori* xc1*g1L*g2R*g0L

! 	
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo

! 	
	  vcAL = vcAL/(4d0*Pi)**2
 	  vcAR = vcAR/(4d0*Pi)**2
 	
 	  vcBL = vcBL/(4d0*Pi)**2
 	  vcBR = vcBR/(4d0*Pi)**2
 	
 	  vcCL = vcCL/(4d0*Pi)**2
 	  vcCR = vcCR/(4d0*Pi)**2
 	
	  vcDL = vcDL/(4d0*Pi)**2
 	  vcDR = vcDR/(4d0*Pi)**2

 	  
c****************************************************
c************ Amplitudenquadrate***********************
c****************************************************

c-----------------------U-channel Neutralino squared------------------------------------------------------- 
	  do neuc=1,6
	  
      n=neuc
      k=neu
      
      colorf = Nc**2
	
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3Lc*g4Rc*
     -       (-(g4R*(muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAL) + 
     -         g4L*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu3*vcAR + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcBL) + g4R*muu3*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) + 
     -      g3Rc*g4Lc*(-(g4L*
     -            (muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAR) + 
     -         g4L*muu3*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -         g4R*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu3*vcAL + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcBR)) + muuqu(k)*
     -       (g3Rc*g4Lc*(-(g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAL)
     -              + g4R*muu4*
     -             (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                  vcAR) + 
     -               4*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g4L*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBR) + 
     -         g3Lc*g4Rc*(-(g4R*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR)
     -              + 2*g4R*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBL + 
     -            g4L*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAL) + 
     -               4*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         sqrtS*(g3Rc*g4Lc*
     -             (g4L*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL
     -                 + g4R*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR
     -               ) + g3Lc*g4Rc*
     -             (g4L*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL
     -                 + g4R*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR
     -               ))*muuqu(k)) + 
     -      muupu(n)*(g3Rc*g4Rc*
     -          (-(g4R*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                 tred)*vcAL) + 
     -            g4L*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            g4R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g3Lc*g4Lc*(-(g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR)
     -              + g4L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            g4R*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAL) + 
     -               2*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         muuqu(k)*(g3Lc*g4Lc*
     -             (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -               2*g4R*muu4*
     -                (-(muu3*vcAR) + 
     -                  sqrtS*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                   vcBL) + 
     -               2*g4L*muu3*sqrtS*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR
     -               ) + g3Rc*g4Rc*
     -             (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -               2*g4R*muu3*sqrtS*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL
     -                 + 2*g4L*muu4*
     -                (-(muu3*vcAL) + 
     -                  sqrtS*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                   vcBR)) + 
     -            sqrtS*(g3Lc*g4Lc*
     -                (g4L*(-1 + muu3**2 + muu4**2)*vcBL + 
     -                  2*g4R*muu3*muu4*vcBR) + 
     -               g3Rc*g4Rc*
     -                (2*g4L*muu3*muu4*vcBL + 
     -                  g4R*(-1 + muu3**2 + muu4**2)*vcBR))*muuqu(k)))
     -      ))/(dconjg(DenGauU(n))*DenGauU(k)))
       
c-----------------------U-channel correction Neutralino with t-channel neutralino------------------------------------------------------- 

      colorf = Nc
	
	  g3Lc = dconjg(ALt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	  g3L = ALt1(neu,isf2,itt2,ftt1,igen2,fgen1)
	
	  g3Rc = dconjg(ARt1(neuc,isf1,itt1,ftt1,igen1,fgen1))
	  g3R = ARt1(neu,isf2,itt2,ftt1,igen2,fgen1)
	
	  g4Lc = dconjg(ALt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  g4L = ALt1(neu,isf1,itt1,ftt2,igen1,fgen2)
	
	  g4Rc = dconjg(ARt1(neuc,isf2,itt2,ftt2,igen2,fgen2))
	  g4R = ARt1(neu,isf1,itt1,ftt2,igen1,fgen2)
	   	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3Lc*g4Rc*
     -       (g4R*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAL + 
     -         g4L*muu4*(muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -             vcAR + sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(-muu1**2 + muu3**2 + tred)*vcBL) - 
     -         g4R*muu3*sqrtS*(muu2**2 - muu4**2 - tred)*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g3Rc*g4Lc*
     -       (g4L*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAR - 
     -         g4L*muu3*sqrtS*(muu2**2 - muu4**2 - tred)*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g4R*muu4*
     -          (muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(-muu1**2 + muu3**2 + tred)*vcBR)) + 
     -      muuqu(neu)*(g3Rc*g4Lc*
     -          (g4L*muu3*(muu2**2 - muu4**2 - tred)*vcAL + 
     -            g4R*muu4*(-((muu3**2 + tred)*vcAR) + 
     -               2*muu3*(-muu2**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBL + muu1**2*(vcAR - 2*muu3*sqrtS*vcBL)) + 
     -            2*g4L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g3Lc*g4Rc*(g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -            2*g4R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g4L*muu4*(-((muu3**2 + tred)*vcAL) + 
     -               2*muu3*(-muu2**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBR + muu1**2*(vcAL - 2*muu3*sqrtS*vcBR))) + 
     -         sqrtS*(g3Rc*g4Lc*
     -             (g4L*muu3*(-muu2**2 + muu4**2 + tred)*vcBL + 
     -               g4R*muu4*(-muu1**2 + muu3**2 + tred)*vcBR) + 
     -            g3Lc*g4Rc*(g4L*muu4*(-muu1**2 + muu3**2 + tred)*
     -                vcBL + g4R*muu3*(-muu2**2 + muu4**2 + tred)*vcBR
     -               ))*muuqu(neu)) + 
     -      muupt(neuc)*(g3Rc*g4Rc*
     -          (-(g4R*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                 tred)*vcAL) + 
     -            g4L*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            g4R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g3Lc*g4Lc*(-(g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR)
     -              + g4L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            g4R*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAL) + 
     -               2*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         muuqu(neu)*(g3Lc*g4Lc*
     -             (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -               2*g4R*muu4*
     -                (-(muu3*vcAR) + 
     -                  sqrtS*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                   vcBL) + 
     -               2*g4L*muu3*sqrtS*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR
     -               ) + g3Rc*g4Rc*
     -             (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -               2*g4R*muu3*sqrtS*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL
     -                 + 2*g4L*muu4*
     -                (-(muu3*vcAL) + 
     -                  sqrtS*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                   vcBR)) + 
     -            sqrtS*(g3Lc*g4Lc*
     -                (g4L*(-1 + muu3**2 + muu4**2)*vcBL + 
     -                  2*g4R*muu3*muu4*vcBR) + 
     -               g3Rc*g4Rc*
     -                (2*g4L*muu3*muu4*vcBL + 
     -                  g4R*(-1 + muu3**2 + muu4**2)*vcBR))*muuqu(neu)))
     -      ))/(dconjg(DenGauT(neuc))*DenGauU(neu)))
      
       
      enddo !Index neuc
       
c-----------------------u-channel correction Neutralino with u-channel gluino------------------------------------------------------- 

      colorf = 0d0
	  
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	  
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  
	  g3L = ALt1(n,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(n,isf2,itt2,ftt1,igen2,fgen1)
	
	  g4L = ALt1(n,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(n,isf1,itt1,ftt2,igen1,fgen2)
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (-(g4R*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -            vcAL) + g4L*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAR) + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL) + 
     -         g4R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g1Lc*g2Rc*
     -       (-(g4R*(muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAL) + 
     -         g4L*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu3*vcAR + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcBL) + g4R*muu3*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) + 
     -      g1Lc*g2Lc*muugl*(-(g4L*muu3*
     -            (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR) + 
     -         g4L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g4R*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR)) + 
     -      g1Rc*g2Lc*(-(g4L*
     -            (muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAR) + 
     -         g4L*muu3*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -         g4R*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu3*vcAL + 
     -            sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -             vcBR)) + muuqu(k)*
     -       (g1Lc*g2Lc*muugl*
     -          (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -            2*g4R*muu4*
     -             (-(muu3*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g4L*muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) +
     -          g1Rc*g2Lc*(-(g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAL)
     -              + g4R*muu4*
     -             (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -                  vcAR) + 
     -               4*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g4L*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBR) + 
     -         g1Lc*g2Rc*(-(g4R*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAR)
     -              + 2*g4R*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBL + 
     -            g4L*muu4*(-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                    tred)*vcAL) + 
     -               4*muu3*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -            2*g4R*muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            2*g4L*muu4*
     -             (-(muu3*vcAL) + 
     -               sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         sqrtS*(g4L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -               g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Lc*g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g1Rc*g2Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + g4R*
     -             (2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -               g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Rc*g2Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g1Lc*g2Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR)*muuqu(k))))/(dconjg(DenGlU)*DenGauU(k)))

     
     
c-----------------------U-channel correction Neutralino with t-channel gluino------------------------------------------------------- 

      colorf = CF*Nc
	  
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)

	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (-(g4R*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -            vcAL) + g4L*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAR) + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL) + 
     -         g4R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g1Lc*g2Rc*
     -       (g4R*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAL + 
     -         g4L*muu4*(muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*
     -             vcAR + sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(-muu1**2 + muu3**2 + tred)*vcBL) - 
     -         g4R*muu3*sqrtS*(muu2**2 - muu4**2 - tred)*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g1Lc*g2Lc*muugl*
     -       (-(g4L*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -            vcAR) + g4L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g4R*muu4*
     -          (-((-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -               vcAL) + 
     -            2*muu3*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR)) + 
     -      g1Rc*g2Lc*(g4L*(muu1**2*muu2**2 - muu3**2*muu4**2 + 
     -            tred - (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             tred + tred**2)*vcAR - 
     -         g4L*muu3*sqrtS*(muu2**2 - muu4**2 - tred)*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g4R*muu4*
     -          (muu3*(muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL + 
     -            sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(-muu1**2 + muu3**2 + tred)*vcBR)) + 
     -      muuqu(k)*(g1Lc*g2Lc*muugl*
     -          (-(g4L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -            2*g4R*muu4*
     -             (-(muu3*vcAR) + 
     -               sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g4L*muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR) +
     -          g1Rc*g2Lc*(g4L*muu3*(muu2**2 - muu4**2 - tred)*
     -             vcAL + g4R*muu4*
     -             (-((muu3**2 + tred)*vcAR) + 
     -               2*muu3*(-muu2**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBL + muu1**2*(vcAR - 2*muu3*sqrtS*vcBL)) + 
     -            2*g4L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(g4R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -            2*g4R*muu3*sqrtS*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL + 
     -            2*g4L*muu4*
     -             (-(muu3*vcAL) + 
     -               sqrtS*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         g1Lc*g2Rc*(g4R*muu3*(muu2**2 - muu4**2 - tred)*vcAR + 
     -            2*g4R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g4L*muu4*(-((muu3**2 + tred)*vcAL) + 
     -               2*muu3*(-muu2**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBR + muu1**2*(vcAL - 2*muu3*sqrtS*vcBR))) + 
     -         sqrtS*(g4L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -               g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Lc*g2Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Rc*g2Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -             g4R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -               g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Rc*g2Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Lc*g2Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)
     -           *muuqu(k))))/(dconjg(DenGlT)*DenGauU(k)))      

      enddo !Index neu

      
c***********************************************************************************************************************************************************
c************ Loop external squark (ifs2,itt2,igen2) & quark (ftt1,fgen1) with Neutralino in u-channel lower vertex ****************************************
c***********************************************************************************************************************************************************	

	  do neu=1,6
	
	  vcAR = 0d0
	  vcAL = 0d0
 	
	  vcBL = 0d0
	  vcBR = 0d0

	  vcCL = 0d0
	  vcCR = 0d0
	
	  vcDL = 0d0
	  vcDR = 0d0
	
c-------------------------------Vertex1 ( Gluon , fermion (l,n), Sfermion (i,k,m) (M0,M1,M2))
	do i=1,2
	   do k=3,4
	     do l=3,4
	      do m=1,3
            do n=1,3
! 
	g0L = ALt1(neu,i,k,l,m,n)
	g0R = ARt1(neu,i,k,l,m,n)
	g1 = BR(1,ftt2,l,fgen2,n)
	g2 = VSfSf(1,isf1,i,itt1,k,igen1,m)
! 	
	Ml0 = 0d0
	Ml1 = Mf(l,n)
	Ml2 = MSf(i,k,m)
! 	
	call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	call myBget(uman,Ml1**2,Ml2**2)	
	
	colori = CF

	vcAL = vcAL + colori*g1*g2*(g0L*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-uman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0R*(muu4*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))
     

	vcAR = vcAR + colori*g1*g2*(g0R*(xb0 + Ml0**2*xc0 + 2*muu1**2*sqrtS**2*xc2 + 
     -     muu4**2*sqrtS**2*xc1+(-uman+muu1**2*sqrtS**2)*(2*xc1+xc2+2*xc0)+muu4*sqrtS*Ml1*(xc2+xc1+2*xc0)) + 
     -     g0L*(muu4*sqrtS*MGau(neu)*(2*xc1+xc2+2*xc0)-MGau(neu)*Ml1*(xc2+2*xc0)))
!      
!      
	vcBL = vcBL + colori*g1*g2*g0R*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	vcBR = vcBR + colori*g1*g2*g0L*(muu4*sqrtS*(2*xc1+xc2+2*xc0)-Ml1*(xc2+2*xc0))
	
	vcCL = vcCL + colori*g1*g2*(g0R*MGau(neu)*(2*xc1+xc2+2*xc0)+g0L*Ml1*(2*xc0+xc1+xc2))
	
	vcCR = vcCR + colori*g1*g2*(g0L*MGau(neu)*(2*xc1+xc2+2*xc0)+g0R*Ml1*(2*xc0+xc1+xc2))
	
	vcDL = vcDL + colori*g1*g2*g0R*(2*xc1+xc2+2*xc0)
	
	vcDR = vcDR + colori*g1*g2*g0L*(2*xc1+xc2+2*xc0)

		enddo
	      enddo
	     enddo
	  enddo
	enddo

c-------------------------------Vertex2 ( Gluino , Sfermion (j,l,n), fermion (k,m) (M0,M1,M2))	
	  do j=1,2
	   do k=3,4
	    do l=3,4
         do m=1,3
		  do n=1,3


	g0L = dconjg(ARt1(neu,j,l,k,n,m))
	g0R = dconjg(ALt1(neu,j,l,k,n,m))
	g1L = GluR(j,l,ftt2,n,fgen2)
	g1R = GluL(j,l,ftt2,n,fgen2)
	g2L = GluR(isf1,itt1,k,igen1,m)
	g2R = GluL(isf1,itt1,k,igen1,m)

	
	Ml0 = MGl
	Ml1 = MSf(j,l,n)
	Ml2 = Mf(k,m)
	
	call myCget(muu4**2*sqrtS**2,uman,muu1**2*sqrtS**2,Ml0**2,Ml1**2,Ml2**2)
	call myBget(uman,Ml1**2,Ml2**2)
	
	colori = CF
! 	
 	vcAL = vcAL -colori*((xb0+Ml0**2*xc0 + muu4**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0L*g2R*g1L
     -	                    + Ml0*Ml2*g0L*g1L*g2L*xc0 - MGau(neu)*muu4*sqrtS*xc1*g0R*g1R*g2L
     -                - Ml2*(muu4*sqrtS*(xc1+xc2)*g1R*g2L*g0L - MGau(neu)*xc2*g1L*g2R*g0R)
     -                - Ml0*(muu4*sqrtS*(xc0+xc1+xc2)*g0L*g1R*g2R - MGau(neu)*(xc2+xc0)*g1L*g2L*g0R))
     
	vcAR = vcAR -colori*((xb0+Ml0**2*xc0 + muu4**2*sqrtS**2*xc1 + muu1**2*sqrtS**2*xc2)*g0R*g2L*g1R
     -	                    + Ml0*Ml2*g0R*g1R*g2R*xc0 - MGau(neu)*muu4*sqrtS*xc1*g0L*g1L*g2R
     -                - Ml2*(muu4*sqrtS*(xc1+xc2)*g1L*g2R*g0R - MGau(neu)*xc2*g1R*g2L*g0L)
     -                - Ml0*(muu4*sqrtS*(xc0+xc1+xc2)*g0R*g1L*g2L - MGau(neu)*(xc2+xc0)*g1R*g2R*g0L))

	vcBL = vcBL - colori* (muu4*sqrtS*xc1*g1R*g2L*g0R+Ml2*xc2*g1L*g2R*g0R+Ml0*(xc2+xc0)*g1L*g2L*g0R)
!      
	vcBR = vcBR - colori* (muu4*sqrtS*xc1*g1L*g2R*g0L+Ml2*xc2*g1R*g2L*g0L+Ml0*(xc2+xc0)*g1R*g2R*g0L)
	
	
	vcCL = vcCL - colori* (MGau(neu)*xc1*g1R*g2L*g0R-ML2*(xc1+xc2)*g1R*g2L*g0L-Ml0*(xc1+xc2+xc0)*g1R*g2R*g0L)
	
	vcCR = vcCL - colori* (MGau(neu)*xc1*g1L*g2R*g0L-ML2*(xc1+xc2)*g1L*g2R*g0R-Ml0*(xc1+xc2+xc0)*g1L*g2L*g0R)
	
	
	vcDL = vcDL + colori* xc1*g1R*g2L*g0R      
	
	vcDR = vcDR + colori* xc1*g1L*g2R*g0L

		  enddo
	     enddo
	    enddo
	   enddo
	  enddo
	
	vcAL = vcAL/(4d0*Pi)**2
 	vcAR = vcAR/(4d0*Pi)**2
 	
 	vcBL = vcBL/(4d0*Pi)**2
 	vcBR = vcBR/(4d0*Pi)**2
 	
 	vcCL = vcCL/(4d0*Pi)**2
 	vcCR = vcCR/(4d0*Pi)**2
 	
	vcDL = vcDL/(4d0*Pi)**2
 	vcDR = vcDR/(4d0*Pi)**2

 	
 	
c****************************************************
c************ Amplitudenquadrate***********************
c****************************************************
 
c-----------------------U-channel Neutralino squared------------------------------------------------------- 
	  do neuc=1,6
	  
      n=neuc
      k=neu
      
      colorf = Nc**2
	
	  g3L = ALt1(neu,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(neu,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(neu,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(neu,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	  
	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3L*g3Lc*g4Rc*
     -       (-((muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAR) + 
     -         muu4*sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -            muu4**2 - tred)*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL) + 
     -      g3L*g3Rc*g4Lc*muu3*
     -       (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -       (-2*muu4*vcAR + 
     -         sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBL)
     -        + g3R*g3Rc*g4Lc*
     -       (-((muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAL) + 
     -         muu4*sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -            muu4**2 - tred)*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR) + 
     -      g3Lc*g3R*g4Rc*muu3*
     -       (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -       (-2*muu4*vcAL + 
     -         sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcBR)
     -        + muuqu(neu)*(-(g3L*
     -            (g3Lc*g4Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3Rc*g4Lc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcAL
     -            ) - g3R*(g3Rc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Lc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*vcAR +
     -          2*g3R*sqrtS*
     -          (2*g3Lc*g4Rc*muu3*muu4*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g3Rc*g4Lc*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred))*vcBL + 
     -         2*g3L*sqrtS*(2*g3Rc*g4Lc*muu3*muu4*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g3Lc*g4Rc*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred))*vcBR + 
     -         sqrtS*(g3L*(g3Lc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Rc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + g3R*
     -             (g3Rc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Lc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR)*muuqu(neu)) + 
     -      muupu(neuc)*(g3Rc*g4Rc*
     -          (-(g3R*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                 tred)*vcAL) + 
     -            g3L*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            g3R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g3Lc*g4Lc*(-(g3L*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR)
     -              + g3L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            g3R*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAL) + 
     -               2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         muuqu(neu)*(-(g3L*
     -               (2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcAL) - 
     -            g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcAR + 
     -            2*g3R*sqrtS*
     -             (g3Rc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Lc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + 2*g3L*sqrtS*
     -             (g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR + sqrtS*
     -             (g3L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcBL + 
     -               g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcBR)*
     -             muuqu(neu)))))/(dconjg(DenGauU(neuc))*DenGauU(neu)))
      
c-----------------------U-t-channel Neutralino squared------------------------------------------------------- 
      
      
      colorf = Nc
		  
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
! 	
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  
! 	  
	  ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g3R*(g3Lc*g4Rc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Rc*g4Lc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2))*vcAL + 
     -      g3L*(g3Rc*g4Lc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Lc*g4Rc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2))*vcAR - 
     -      g3L*sqrtS*(g3Lc*g4Rc*muu4*(muu1**2 - muu3**2 - tred) + 
     -         g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -       (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*vcBL -
     -        g3R*sqrtS*(g3Rc*g4Lc*muu4*(muu1**2 - muu3**2 - tred) + 
     -         g3Lc*g4Rc*muu3*(muu2**2 - muu4**2 - tred))*
     -       (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*vcBR +
     -        muuqu(k)*(g3L*
     -          (g3Lc*g4Rc*muu4*(muu1**2 - muu3**2 - tred) + 
     -            g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*vcAL + 
     -         g3R*g3Rc*g4Lc*
     -          (muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            2*sqrtS*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL) + 
     -         g3Lc*g3R*g4Rc*muu3*
     -          (-((muu4**2 + tred)*vcAR) + 
     -            2*muu4*(-muu1**2 + muu3**2 + muu4**2)*sqrtS*vcBL + 
     -            muu2**2*(vcAR - 2*muu4*sqrtS*vcBL)) + 
     -         2*g3L*sqrtS*(g3Rc*g4Lc*muu3*muu4*
     -             (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g3Lc*g4Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))*vcBR + 
     -         sqrtS*(g3L*(g3Lc*g4Rc*muu4*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -             g3R*(g3Rc*g4Lc*muu4*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)
     -           *muuqu(k)) + 
     -      muupt(n)*(g3Rc*g4Rc*
     -          (-(g3R*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                 tred)*vcAL) + 
     -            g3L*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            g3R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g3Lc*g4Lc*(-(g3L*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR)
     -              + g3L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL + 
     -            g3R*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAL) + 
     -               2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         muuqu(k)*(-(g3L*
     -               (2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcAL) - 
     -            g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcAR + 
     -            2*g3R*sqrtS*
     -             (g3Rc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Lc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + 2*g3L*sqrtS*
     -             (g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR + sqrtS*
     -             (g3L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*vcBL + 
     -               g3R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*vcBR)*
     -             muuqu(k)))))/(dconjg(DenGauT(n))*DenGauU(k)))
       
      enddo !End of neuc
       
c-----------------------U-Neutralino t-channel Gluino squared------------------------------------------------------- 
      
      
      colorf  = CF*Nc      
      
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
      
      
      ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal((sqrtS**4*(g1Lc*(g2Rc*g3R*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAL - 
     -         g2Lc*g3R*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*vcAL - 
     -         g2Lc*g3L*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR + 
     -         g2Rc*g3L*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)*vcAR + 
     -         g2Lc*g3L*(-1 + muu3**2 + muu4**2)*muugl*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g2Rc*g3L*muu4*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-muu1**2 + muu3**2 + tred)*vcBL + 
     -         g3R*muu3*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-muu2**2 + muu4**2 + tred))*vcBR) + 
     -      g1Rc*(g2Rc*muugl*
     -          (-(g3R*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -                 tred)*vcAL) + 
     -            g3L*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAR) + 
     -               2*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            g3R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR) + 
     -         g2Lc*(g3R*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -               tred**2)*vcAL + 
     -            g3L*muu3*(muu4*
     -                (muu1**2 + muu2**2 - muu3**2 - muu4**2)*vcAR + 
     -               sqrtS*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred)*(-muu2**2 + muu4**2 + tred)*
     -                vcBL) - 
     -            g3R*muu4*sqrtS*(muu1**2 - muu3**2 - tred)*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR)) + 
     -      muuqu(k)*(g1Lc*g2Lc*muugl*
     -          (-(g3L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -            2*g3R*muu3*
     -             (-(muu4*vcAR) + 
     -               sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g3L*muu4*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR) +
     -          g1Lc*g2Rc*(g3L*muu4*(muu1**2 - muu3**2 - tred)*
     -             vcAL + g3R*muu3*
     -             (-((muu4**2 + tred)*vcAR) + 
     -               2*muu4*(-muu1**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBL + muu2**2*(vcAR - 2*muu4*sqrtS*vcBL)) + 
     -            2*g3L*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBR) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(g3R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -            2*g3R*muu4*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*g3L*muu3*
     -             (-(muu4*vcAL) + 
     -               sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcBR)) + 
     -         g1Rc*g2Lc*(g3R*muu4*(muu1**2 - muu3**2 - tred)*vcAR + 
     -            2*g3R*sqrtS*
     -             ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))*vcBL + 
     -            g3L*muu3*(-((muu4**2 + tred)*vcAL) + 
     -               2*muu4*(-muu1**2 + muu3**2 + muu4**2)*sqrtS*
     -                vcBR + muu2**2*(vcAL - 2*muu4*sqrtS*vcBR))) + 
     -         sqrtS*(g3L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -               g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Lc*g2Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Rc*g2Lc*muu3*(-muu2**2 + muu4**2 + tred))*vcBL +
     -             g3R*(2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -               g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Rc*g2Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g1Lc*g2Rc*muu3*(-muu2**2 + muu4**2 + tred))*vcBR)
     -           *muuqu(k))))/(dconjg(DenGlT)*DenGauU(k)))


c-----------------------U-Neutralino u-channel Gluino squared------------------------------------------------------- 

      colorf = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(k,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(k,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(k,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(k,isf1,itt1,ftt2,igen1,fgen2))
	  
	  
      ststQQ_M2VertexT = ststQQ_M2VertexT + 2d0*colorf*dreal( (sqrtS**4*(g1Rc*g2Rc*muugl*
     -       (-(g3R*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*
     -            vcAL) + g3L*muu3*
     -          (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAR) + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBL) + 
     -         g3R*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBR) + g1Rc*g2Lc*
     -       (-(g3R*(muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAL) + 
     -         g3L*muu3*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu4*vcAR + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcBL) + g3R*muu4*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR) + 
     -      g1Lc*g2Lc*muugl*(-(g3L*muu4*
     -            (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR) + 
     -         g3L*(-1 + muu3**2 + muu4**2)*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          vcBL + g3R*muu3*
     -          (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -               vcAL) + 
     -            2*muu4*sqrtS*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*vcBR)) + 
     -      g1Lc*g2Rc*(-(g3L*
     -            (muu3**4 - muu4**2 + 
     -              muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -              (muu4**2 - tred)**2 + 
     -              muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -              muu2**2*tred)*vcAR) + 
     -         g3L*muu4*sqrtS*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -         g3R*muu3*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred)*(-2*muu4*vcAL + 
     -            sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -             vcBR)) + muuqu(k)*
     -       (g1Lc*g2Lc*muugl*
     -          (-(g3L*(-1 + muu3**2 + muu4**2)*vcAL) + 
     -            2*g3R*muu3*
     -             (-(muu4*vcAR) + 
     -               sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g3L*muu4*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBR) +
     -           g1Lc*g2Rc*(-(g3L*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAL)
     -              + g3R*muu3*
     -             (-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -                  vcAR) + 
     -               4*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBL) + 
     -            2*g3L*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBR) + 
     -         g1Rc*g2Lc*(-(g3R*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcAR)
     -              + 2*g3R*sqrtS*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)*vcBL + 
     -            g3L*muu3*(-((-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)*vcAL) + 
     -               4*muu4*sqrtS*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*vcBR)) + 
     -         g1Rc*g2Rc*muugl*
     -          (-(g3R*(-1 + muu3**2 + muu4**2)*vcAR) + 
     -            2*g3R*muu4*sqrtS*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)*vcBL + 
     -            2*g3L*muu3*
     -             (-(muu4*vcAL) + 
     -               sqrtS*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                  tred)*vcBR)) + 
     -         sqrtS*(g3L*(2*g1Rc*g2Rc*muu3*muu4*muugl + 
     -               g1Lc*g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Lc*g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g1Rc*g2Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBL + g3R*
     -             (2*g1Lc*g2Lc*muu3*muu4*muugl + 
     -               g1Rc*g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g1Rc*g2Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g1Lc*g2Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             vcBR)*muuqu(k))))/(dconjg(DenGlU)*DenGauU(k)))
     
     
      enddo !Enddo of neu


	end
