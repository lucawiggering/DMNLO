* All loop functions are all defined with a prefactor c_eps = Gamma(1+eps) (4pi)^eps
* which when multiplied with the single UV pole give the minimal subtraction Delta_UV
*
* With c_eps factor - UVdiv stands for the pole only (not for Delta_UV anymore)
*                   - IRdiv stands also for the pole only (pole in the same eps as UVdiv)
* eps = (4-D)/2 

*************************************************************
	subroutine myAget(m02)
*************************************************************

#include "../DMNLO_LoopIntegrals.h"

	double precision m02
* functions
	double precision A01

	xa0 = A01(dsqrt(m02))

	end

*************************************************************
	subroutine myBget(p1,m02,m12)
*************************************************************

#include "../DMNLO_LoopIntegrals.h"

	double precision m02,m12
	double precision p1,m0,m1
	double precision f1D
* functions
	complex*16 B02
	double precision A01

c UV single pole
      if (EpsPole.eq.1d0) then
!--------------------------------------------------------------------            
       xb0  = 1d0*UVdiv
       xb1  = -0.5d0*UVdiv
       xb00 = -1d0*(p1/12d0-(m02+m12)/4d0)*UVdiv
       xb11 = 1d0/3d0*UVdiv

c IR single pole
      else if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------            
       
       if ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.eq.0d0)) then

	   xb0  = B02(p1,dsqrt(m02),dsqrt(m12))
	   xb1  = -1/2d0*B02(p1,dsqrt(m02),dsqrt(m12))
 	   xb00 = 0d0      
 	   xb11 = 1/3d0*B02(p1,dsqrt(m02),dsqrt(m12))

       else 
         xb0  = 0d0
         xb1  = 0d0
		 xb00 = 0d0      
		 xb11 = 0d0

       endif

c regular part
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------                  
	xb0 = B02(p1,dsqrt(m02),dsqrt(m12))

*********************************
* B^{mu} coefficients 
*
	f1D = p1 - m12 + m02

**** cases for special arguments

c b1(p1,m0,m1)
	if (p1.ne.0d0) then
	  xb1 = (1d0/(2d0*p1))*(A01(dsqrt(m02))-A01(dsqrt(m12))-f1D*xb0) 	
c b1(0,m0,m1)
	elseif ((p1.eq.0d0).and.(m02.ne.m12)) then
	  xb1 = A01(dsqrt(m12))/(m02 - m12) + (-m02**2/4d0 + m12**2/4d0 - (m02*A01(dsqrt(m02)))/2d0 + (m12*A01(dsqrt(m12)))/2d0)/(m02 - m12)**2
c b1(0,m0,m0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.ne.0d0)) then
	  xb1 = -1/2d0*B02(p1,dsqrt(m02),dsqrt(m12))
c b1(0,0,0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.eq.0d0)) then
	  xb1 = -1/2d0*B02(p1,dsqrt(m02),dsqrt(m12))

	endif

*********************************
* B^{mu nu} coefficients 
*

	xb00 = (m02 + 2*m02*xb0 + m12 - p1/3d0 + xb1*(m02 - m12 + p1) + A01(dsqrt(m12)))/6d0

**** cases for special arguments

c b11(p1,m0,m1)
	if (p1.ne.0d0) then
	  xb11 = (-m02 - 2d0*m02*xb0 - m12 + p1/3d0 - 4*xb1*(m02 - m12 + p1) + 2*A01(dsqrt(m12)))/(6d0*p1)
c b11(0,m0,m1)
	elseif ((p1.eq.0d0).and.(m02.ne.m12)) then
	  xb11 = -(A01(dsqrt(m12))/(m02 - m12)) - (2*m02*
     -     (-m02**2/4d0 + m12**2/4d0 - (m02*A01(dsqrt(m02)))/2d0 + (m12*A01(dsqrt(m12)))/2d0))/(m02 - m12)**3 + 
     -  (2*(-m02**3/9d0 + m12**3/9d0 - (m02**2*A01(dsqrt(m02)))/3d0 + (m12**2*A01(dsqrt(m12)))/3d0))/(m02 - m12)**3
c b11(0,m0,m0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.ne.0d0)) then
	  xb11 = 1/3d0*B02(p1,dsqrt(m02),dsqrt(m12))
c b11(0,0,0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.eq.0d0)) then
	  xb00 = 0d0      
	  xb11 = 1/3d0*B02(p1,dsqrt(m02),dsqrt(m12))

	endif

!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        xb0  = 0d0
        xb1  = 0d0
        xb00 = 0d0
        xb11 = 0d0
      endif
      

	end

*************************************************************
	subroutine myBPget(p1,m02,m12)
*************************************************************

#include "../DMNLO_LoopIntegrals.h"

	double precision m02,m12
	double precision p1,m0,m1
	double precision f1D
* swap
	complex*16 swb0,swb1,swb00,swb11
	complex*16 pb0,pb1,pb00,pb11
* functions
	complex*16 B02,B02p
	double precision A01

c IR single pole
      if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------            
       if (((p1.eq.m02).and.(m02.ne.0.d0).and.(m12.eq.0.d0)).or.
     &   ((p1.eq.m12).and.(m12.ne.0.d0).and.(m02.eq.0.d0)) ) then

        pb0  = 0d0
        pb1  = 0d0
        pb00 = 0d0
        pb11 = 0d0
        
        f1D = p1 - m12 + m02
  	
	  xbp0 = B02p(p1,dsqrt(m02),dsqrt(m12))
	  xbp1 = -(1d0/p1)*(pb1 + 1/2d0*pb0 + 1/2d0*f1D*xbp0) 	
	  xbp00 = (2d0*m02*xbp0 + xb1 + xbp1*(m02 - m12 + p1))/6d0
  	  xbp11 = (-6d0*pb11 - 2d0*m02*xbp0 - 4*xbp1*(m02 - m12 + p1) - 4*xb1)/(6d0*p1)
	
       else

		xbp0  = 0d0
		xbp1  = 0d0
        xbp00 = 0d0
        xbp11 = 0d0
        
       endif
      
c regular part
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------                  
	xbp0 = B02p(p1,dsqrt(m02),dsqrt(m12))

*********************************
* Bp^{mu} coefficients 
*

	swb0  = xb0
	swb1  = xb1
	swb00 = xb00
	swb11 = xb11

	call myBget(p1,m02,m12)
	pb0  = xb0
	pb1  = xb1
	pb00 = xb00
	pb11 = xb11

	xb0  = swb0
	xb1  = swb1
	xb00 = swb00
	xb11 = swb11
*
**** cases for special arguments

	f1D = p1 - m12 + m02

c b1'(p1,m0,m1)
	if (p1.ne.0d0) then
	  xbp1 = -(1d0/p1)*(pb1 + 1/2d0*pb0 + 1/2d0*f1D*xbp0) 	
c b1'(0,m0,m1)
	elseif ((p1.eq.0d0).and.(m02.ne.m12).and.(m02.ne.0d0).and.(m12.ne.0d0)) then
	  xbp1 = -1/(6.*(-m02 + m12)) + (m02*(m02 + m12))/(2.*(-m02 + m12)**3) - (m02**2*m12*dlog(m12/m02))/(-m02 + m12)**4
c b1'(0,m0,m0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.ne.0d0)) then
	  xbp1 = -1/(12d0*m02)
c b1'(0,m0,0)
	elseif ((p1.eq.0d0).and.(m02.ne.0d0).and.(m12.eq.0d0)) then
	  xbp1 = -1/(3d0*m02)
c b1'(0,0,m1)
	elseif ((p1.eq.0d0).and.(m12.ne.0d0).and.(m02.eq.0d0)) then
	  xbp1 = -1/(6d0*m12)
c b1'(0,0,0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.eq.0d0)) then
	  xbp1 = 0d0

	endif

*********************************
* Bp^{mu nu} coefficients 
*

	xbp00 = (- 1/3d0 + 2d0*m02*xbp0 + xb1 + xbp1*(m02 - m12 + p1))/6d0

c b11'(p1,m0,m1)
	if (p1.ne.0d0) then
	  xbp11 = (-6d0*pb11 - 2d0*m02*xbp0 + 1/3d0 - 4*xbp1*(m02 - m12 + p1) - 4*xb1)/(6d0*p1)
c b11'(0,m0,m1)
	elseif ((p1.eq.0d0).and.(m02.ne.m12).and.(m02.ne.0d0).and.(m12.ne.0d0)) then
	  xbp11 = ((m02 - m12)*(3*m02**3 + 13*m02**2*m12 - 5*m02*m12**2 + m12**3) + 12*m02**3*m12*dlog(m12/m02))/
     -  (12.*(m02 - m12)**5)
c b11'(0,m0,m0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.ne.0d0)) then
	  xbp11 = 1/(20d0*m02)
c b11'(0,m0,0)
	elseif ((p1.eq.0d0).and.(m02.ne.0d0).and.(m12.eq.0d0)) then
	  xbp11 = 1/(4d0*m02)
c b11'(0,0,m1)
	elseif ((p1.eq.0d0).and.(m12.ne.0d0).and.(m02.eq.0d0)) then
	  xbp11 = 1/(12d0*m12)
c b11'(0,0,0)
	elseif ((p1.eq.0d0).and.(m02.eq.m12).and.(m02.eq.0d0)) then
	  xbp00 = -1/12d0*(pb0-1d0)
	  xbp11 = 0d0

	endif
!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        xbp0  = 0d0
        xbp1  = 0d0
        xbp00 = 0d0
        xbp11 = 0d0
      endif

	end

*************************************************************
	subroutine myCget(p1,p21,p2,m02,m12,m22)
*************************************************************
      implicit none
#include "../DMNLO_LoopIntegrals.h"

	double precision m02,m12,m22
	double precision p1,p2,p21,m0,m1,m2
	double precision X2(2,2), detX2
	double precision f1D, f2D, test
	complex*16 R1D, R2D, R3D, R4D, R5D, R6D, sC03, x0
* swap
	complex*16 swb0,swb1,swb00,swb11
	complex*16 cb0_1,cb1_1,cb00_1,cb11_1,cb0_2,cb1_2,cb00_2,cb11_2,cb0_3,cb1_3,cb00_3,cb11_3
* functions
	complex*16 B02,C03f
	double precision A01
* subroutines
	external myBget, alternativeCget

!	ieps = dcmplx(0.D0,1.D-17)
	
	xc0 = C03f(p1,p21,p2,dsqrt(m02),dsqrt(m12),dsqrt(m22))
	
cc	write(*,*)'Incoming C0',xc0
cc	write(*,*)'C arguments',p1,p21,p2,m02,m12,m22

	
! Special cases for alternative C tensor determination relevant for DD -- Work in progress by PS!
! t-channel with both momenta equal and zero momentum transfer C(p**2, 0, p**2, m0**2, m1**2, m2**2)	
	if ((p21.eq.0d0).and.(p1.eq.p2)) then

!		write(*,*)'WARNING! Encountered t-channel critical for tensor reduction'
!		write(*,*)''
		call alternativeCget(p1,p21,p2,m02,m12,m22)

! s-channel or u-channel with one momentum zero and zero momentum transfer C(p**2, p**2, 0, m0**2, m1**2, m2**2)
	elseif((p21.eq.p1).and.(p2.eq.0d0)) then

!		write(*,*)'WARNING! Encountered s-channel or u-channel critical for tensor reduction'
!		write(*,*)''
		call alternativeCget(p1,p21,p2,m02,m12,m22)

! s-channel or u-channel with one momentum zero and zero momentum transfer C(0, p**2, p**2, m0**2, m1**2, m2**2)
	elseif((p21.eq.p2).and.(p1.eq.0d0)) then

!		write(*,*)'WARNING! Encountered s-channel or u-channel critical for tensor reduction'
!		write(*,*)''
		call alternativeCget(p1,p21,p2,m02,m12,m22)

! s-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)+m(p2))**2, p2**2, m0**2, m1**2, m2**2)
	elseif(p21.eq.((dsqrt(p1) + dsqrt(p2))**2).or.(abs(p21 - (dsqrt(p1) + dsqrt(p2))**2).lt.0.001d0)) then

!		write(*,*)'WARNING! Encountered s-channel critical for tensor reduction'
!		write(*,*)''
		call alternativeCget(p1,p21,p2,m02,m12,m22)

! u-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)-m(p2))**2, p2**2, m0**2, m1**2, m2**2)
	elseif(p21.eq.((dsqrt(p1) - dsqrt(p2))**2).or.(abs(p21 - (dsqrt(p1) - dsqrt(p2))**2).lt.0.001d0)) then

!		write(*,*)'WARNING! Encountered u-channel critical for tensor reduction'
!		write(*,*)''
		call alternativeCget(p1,p21,p2,m02,m12,m22)

! Regular case
	else
	
*********************************
* C^{mu} coefficients 
*

	detX2 = (-p1**2 - (p2 - p21)**2 + 2*p1*(p2 + p21))/4d0

cc	write(*,*)'Gram Determinante = ',detX2
	
	X2(1,1) = p2
	X2(2,2) = p1
	X2(1,2) = (-p1 - p2 + p21)/2d0
	X2(2,1) = X2(1,2)

	f1D = p1 - m12 + m02
	f2D = p2 - m22 + m02

	R1D = 0.5d0*(B02(p2,dsqrt(m02),dsqrt(m22)) - B02(p21,dsqrt(m12),dsqrt(m22)) - f1D*xc0)
	R2D = 0.5d0*(B02(p1,dsqrt(m02),dsqrt(m12)) - B02(p21,dsqrt(m12),dsqrt(m22)) - f2D*xc0)

****
	xc1 = 1d0/detX2*(X2(1,1)*R1D + X2(1,2)*R2D) 	
	xc2 = 1d0/detX2*(X2(2,1)*R1D + X2(2,2)*R2D)

*********************************
* C^{munu} coefficients 
*
******** B coefficients

	swb0  = xb0
	swb1  = xb1
	swb00 = xb00
	swb11 = xb11

	call myBget(p21,m12,m22)
	cb0_1  = xb0
	cb1_1  = xb1
	cb00_1  = xb00
	cb11_1  = xb11
	call myBget(p2,m02,m22)
	cb0_2  = xb0
	cb1_2  = xb1
	cb00_2  = xb00
	cb11_2  = xb11
	call myBget(p1,m02,m12)
	cb0_3  = xb0
	cb1_3  = xb1
	cb00_3  = xb00
	cb11_3  = xb11

	xb0  = swb0
	xb1  = swb1
	xb00 = swb00
	xb11 = swb11
*********

    
c C0s give no additional finite parts
c
c c00(0,s,0;0,0,0) - double pole (soft & collinear)
c c0(0,s,q3^2;0,0,q3) - double pole (soft & collinear)
c c0(0,s,k2^2;0,0,q3) - single pole (pure collinear)
c c00(q2^2,s,q3^2;0,q2,q3) - single pole (pure soft)
c
c therefore treat all similar:
          

c	xc00 = 1/4d0*(cb0_1 + 1 + 2*m02*xc0 + f1D*xc1 + f2D*xc2)	

	xc00 = 1/4d0*(cb0_1 + 2*m02*xc0 + f1D*xc1 + f2D*xc2)	
      if (EpsPole.eq.0d0) then 
            xc00 = xc00 + 1/4d0
      endif
	
	R3D = 0.5d0*(cb0_1 + cb1_1 - f1D*xc1) - xc00
	R4D = 0.5d0*(cb0_1 + cb1_1 + cb1_3 - f2D*xc1)
	R5D = 0.5d0*(cb1_2 - cb1_1 - f1D*xc2)
	R6D = -0.5d0*(cb1_1 + f2D*xc2) - xc00

****
	xc11 = 1d0/detX2*(X2(1,1)*R3D + X2(1,2)*R4D) 	
	xc12 = 1d0/detX2*(X2(2,1)*R3D + X2(2,2)*R4D) 	
	xc22 = 1d0/detX2*(X2(2,1)*R5D + X2(2,2)*R6D) 	

*********************************
* C^{mu nu rho} coefficients 
*
*********************************

	xc001 = 1/6d0*(-cb0_1 - cb1_1 + f2D*xc12 + 2d0*m02*xc1 + f1D*xc11)
	xc002 = 1/6d0*(       + cb1_1 + f1D*xc12 + 2d0*m02*xc2 + f2D*xc22)
	
	if (EpsPole.eq.0d0) then 
	      xc001 = xc001 - 1/18d0
	      xc002 = xc002 - 1/18d0
	endif

	R3D = 0.5d0*(-cb0_1 - 2d0*cb1_1 - cb11_1 - xc11*f1D - 4d0*xc001)
	R4D = 0.5d0*(-cb0_1 - 2d0*cb1_1 - cb11_1 - xc11*f2D + cb11_3)

	R5D = 0.5d0*(-cb11_1 + cb11_2 - xc22*f1D)
	R6D = 0.5d0*(-cb11_1           -xc22*f2D - 4d0*xc002)

	xc111 = 1d0/detX2*(X2(1,1)*R3D + X2(1,2)*R4D) 	
	xc112 = 1d0/detX2*(X2(2,1)*R3D + X2(2,2)*R4D) 	

	xc122 = 1d0/detX2*(X2(1,1)*R5D + X2(1,2)*R6D) 	
	xc222 = 1d0/detX2*(X2(2,1)*R5D + X2(2,2)*R6D) 	

! End if for alternative C tensor determination	
	endif	
	
	end


*************************************************************
	subroutine myDget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32)
*************************************************************

#include "../DMNLO_LoopIntegrals.h"
	
	double precision m02,m12,m22,m32
	double precision p1,p2,p3,p21,p32,p13,m0,m1,m2,m3
	double precision X3(3,3),detX3
	double precision f1D, f2D, f3D
	complex*16 R1D, R2D, R3D, R4D, R5D, R6D, R7D, R8D, R9D, R10D, R11D, R12D
	complex*16 R01D,R02D,R03D,R13D,R21D,R22D,R23D,R31D,R32D,R33D
	complex*16 R51D,R52D,R53D,R61D,R62D,R63D,R91D,R92D,R93D
	double precision accthr, accthr1
	double precision MACHEPS
           
* swap
	complex*16 swc0,swc1,swc2,swc00,swc11,swc12,swc22
	complex*16 dc0_0,dc1_0,dc2_0,dc00_0,dc11_0,dc12_0,dc22_0
	complex*16 dc0_1,dc1_1,dc2_1,dc00_1,dc11_1,dc12_1,dc22_1
	complex*16 dc0_2,dc1_2,dc2_2,dc00_2,dc11_2,dc12_2,dc22_2
	complex*16 dc0_3,dc1_3,dc2_3,dc00_3,dc11_3,dc12_3,dc22_3
	complex*16 tempD1, tempD2, tempD3, tempD00, tempD11, tempD12, tempD13, tempD22, tempD23, tempD33
	
* functions
	complex*16 C03f,D04f
* subroutines
	external myCget, alternativeDget

	xd0 = 0d0
	xd1 = 0d0
	xd2 = 0d0
	xd3 = 0d0

	R1D = 0d0
	R2D = 0d0
	R3D = 0d0
      
***************************************** 
*  C-coefficients needed in D-reduction
*
*****************************************
	swc0  = xc0
	swc1  = xc1
	swc2  = xc2
	swc00 = xc00
	swc11 = xc11
	swc12 = xc12
	swc22 = xc22

	call myCget(p21,p32,p13,m12,m22,m32)
	dc0_0  = xc0
	dc1_0  = xc1
	dc2_0  = xc2
	dc00_0 = xc00
	dc11_0 = xc11
	dc12_0 = xc12
	dc22_0 = xc22
	call myCget(p2,p32,p3,m02,m22,m32)
	dc0_1  = xc0
	dc1_1  = xc1
	dc2_1  = xc2
	dc00_1 = xc00
	dc11_1 = xc11
	dc12_1 = xc12
	dc22_1 = xc22
	call myCget(p1,p13,p3,m02,m12,m32)
	dc0_2  = xc0
	dc1_2  = xc1
	dc2_2  = xc2
	dc00_2 = xc00
	dc11_2 = xc11
	dc12_2 = xc12
	dc22_2 = xc22
	call myCget(p1,p21,p2,m02,m12,m22)
	dc0_3  = xc0
	dc1_3  = xc1
	dc2_3  = xc2
	dc00_3 = xc00
	dc11_3 = xc11
	dc12_3 = xc12
	dc22_3 = xc22

	xc0  = swc0
	xc1  = swc1
	xc2  = swc2
	xc00 = swc00
	xc11 = swc11
	xc12 = swc12
	xc22 = swc22

**************************************
* D_0 basic scalar function

	xd0 = D04f(p1,p21,p32,p3,p2,p13,dsqrt(m02),dsqrt(m12),dsqrt(m22),dsqrt(m32))

*********************************
* D^{mu} coefficients 
*	
*********************************

	detX3 = (-(p13**2*p2) - p1**2*p32 + p1*((p13 - p21)*(p2 - p3) + (p13 + p2 + p21 + p3)*p32 - p32**2) - p3*(p21*(p21 + p3 - p32) + 
     -    p2*(-p21 + p32)) + p13*(-p2**2 + p21*(p3 - p32) + p2*(p21 + p3 + p32)) )/4d0
	
! Setting numerical noise to zero
	if(abs(detX3).lt.0.1) then
!		write(*,*)'WARNING! |Gramdet| < 0.1, automatically set to 0. Old detX3 = ', detX3
		detX3 = 0d0
	endif

! Special cases for alternative D tensor determination relevant for DD -- Work in progress by PS!

! Box B2 or B4 with zero momentum transfer D(p, p, p3, m0**2, m1**2, m2**2, m3**2)
	if ((p1.eq.p2).and.(p21.eq.0d0).and.(p13.eq.p32).and.(p1.ne.p3)) then

!		write(*,*)'WARNING! Encountered Box B2 or B4 critical for D Tensor reduction'
!		write(*,*)''

		call alternativeDget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32)

! Box B1 or B3 with zero momentum transfer D(p1, p, p, m0**2, m1**2, m2**2, m3**2)
	elseif ((p3.eq.p2).and.(p32.eq.0d0).and.(p13.eq.p21).and.(p1.ne.p3)) then

!		write(*,*)'WARNING! Encountered Box B1 or B3 critical for D Tensor reduction'
!		write(*,*)''

!Rotate arguments so that the case becomes equivalent to the previous one (1 <-> 3)
		call alternativeDget(p3,p32,p21,p1,p2,p13,m02,m32,m22,m12)

		D1temp = xd1
		D2temp = xd2
		D3temp = xd3

		D00temp = xd00
		D11temp = xd11
		D12temp = xd12
		D13temp = xd13
		D22temp = xd22
		D23temp = xd23
		D33temp = xd33

		xd1 = D3temp
		xd2 = D2temp
		xd3 = D1temp

		xd00 = D00temp
		xd11 = D33temp
		xd12 = D23temp
		xd13 = D13temp
		xd22 = D22temp
		xd23 = D12temp
		xd33 = D11temp

! Regular case
	else
		
	X3(1,1) = (-p2**2 - (p3 - p32)**2 + 2*p2*(p3 + p32))/4.
	X3(2,2) = (-p1**2 - (p13 - p3)**2 + 2*p1*(p13 + p3))/4.
	X3(3,3) = (-p1**2 - (p2 - p21)**2 + 2*p1*(p2 + p21))/4.

	X3(1,2) = (p1*(p2 - p3 - p32) - p13*(p2 + p3 - p32) + p3*(-p2 + 2*p21 + p3 - p32))/4.
	X3(1,3) = (2*p13*p2 + p2**2 - p2*p21 - p2*p3 - p21*p3 - p2*p32 + p21*p32 - p1*(p2 - p3 + p32))/4.
	X3(2,3) = (p1**2 - (p2 - p21)*(p13 - p3) - p1*(p13 + p2 + p21 + p3 - 2*p32))/4.

	X3(2,1) = X3(1,2)  
	X3(3,1) = X3(1,3)  
	X3(3,2) = X3(2,3)  

	f1D = p1 - m12 + m02
	f2D = p2 - m22 + m02
	f3D = p3 - m32 + m02

	R1D = 0.5d0*(dc0_1 - dc0_0 - f1D*xd0)
	R2D = 0.5d0*(dc0_2 - dc0_0 - f2D*xd0)
	R3D = 0.5d0*(dc0_3 - dc0_0 - f3D*xd0)

c to remove IR cancellations in tensor coefficients
      if (EpsPole.lt.0d0) then
            accthr = dabs(real(1d3*(f1D+f2D+f3D)*xd0*MACHEPS()))
            if (dabs(dreal(R1D)).lt.accthr) R1D = 0d0
            if (dabs(dreal(R2D)).lt.accthr) R2D = 0d0
            if (dabs(dreal(R3D)).lt.accthr) R3D = 0d0
      endif

**** D1, D2, D3 coefficients

	xd1 = 1d0/detX3*(X3(1,1)*R1D + X3(1,2)*R2D + X3(1,3)*R3D) 	
	xd2 = 1d0/detX3*(X3(2,1)*R1D + X3(2,2)*R2D + X3(2,3)*R3D) 	
	xd3 = 1d0/detX3*(X3(3,1)*R1D + X3(3,2)*R2D + X3(3,3)*R3D) 	

c to remove IR cancellations in tensor coefficients
      if (EpsPole.lt.0d0) then
            accthr1 = accthr*(dabs(X3(1,1)) + dabs(X3(1,2)) + dabs(X3(1,3)) +
     &                        dabs(X3(2,1)) + dabs(X3(2,2)) + dabs(X3(2,3)) +
     &                        dabs(X3(3,1)) + dabs(X3(3,2)) + dabs(X3(3,3)))/(3d0*dabs(detX3))
            
            if (dabs(real(xd1)).lt.accthr1) xd1 = 0d0
            if (dabs(real(xd2)).lt.accthr1) xd2 = 0d0
            if (dabs(real(xd3)).lt.accthr1) xd3 = 0d0
      endif

*********************************
*    D^{mu nu} coefficients 
*	
*********************************

**** D00 coefficient ( no finite parts due to IR divergence)

	xd00 = 1/2d0*(dc0_0 + 2*m02*xd0 + f1D*xd1 + f2D*xd2 + f3D*xd3)	

*******
	
	R4D  = 0.5d0*(        dc0_0 + dc1_0 + dc2_0 - f1D*xd1) - xd00
	R5D  = 0.5d0*(dc1_2 + dc0_0 + dc1_0 + dc2_0 - f2D*xd1)
	R6D  = 0.5d0*(dc1_3 + dc0_0 + dc1_0 + dc2_0 - f3D*xd1)
	R7D  = 0.5d0*(dc1_1 - dc1_0                 - f1D*xd2)
	R8D  = 0.5d0*(      - dc1_0                 - f2D*xd2) - xd00
	R9D  = 0.5d0*(dc2_3 - dc1_0                 - f3D*xd2)
	R10D = 0.5d0*(dc2_1 - dc2_0                 - f1D*xd3)
	R11D = 0.5d0*(dc2_2 - dc2_0                 - f2D*xd3)
	R12D = 0.5d0*(      - dc2_0                 - f3D*xd3) - xd00

**** D11, D12, D13, D22, D23, D33 coefficients

	xd11 = 1d0/detX3*(X3(1,1)*R4D + X3(1,2)*R5D + X3(1,3)*R6D) 	
	xd12 = 1d0/detX3*(X3(2,1)*R4D + X3(2,2)*R5D + X3(2,3)*R6D) 	
	xd13 = 1d0/detX3*(X3(3,1)*R4D + X3(3,2)*R5D + X3(3,3)*R6D) 	

c alternative xd12
c	xd12 = 1d0/detX3*(X3(1,1)*R7D + X3(1,2)*R8D + X3(1,3)*R9D) 	
	xd22 = 1d0/detX3*(X3(2,1)*R7D + X3(2,2)*R8D + X3(2,3)*R9D) 	
	xd23 = 1d0/detX3*(X3(3,1)*R7D + X3(3,2)*R8D + X3(3,3)*R9D) 	

c alternative xd13, xd23
c	xd13 = 1d0/detX3*(X3(1,1)*R10D + X3(1,2)*R11D + X3(1,3)*R12D) 	
c	xd23 = 1d0/detX3*(X3(2,1)*R10D + X3(2,2)*R11D + X3(2,3)*R12D) 	
	xd33 = 1d0/detX3*(X3(3,1)*R10D + X3(3,2)*R11D + X3(3,3)*R12D) 	

c to remove IR cancellations in tensor coefficients
      if (EpsPole.lt.0d0) then
            accthr = dabs(1d3*(f1D+f2D+f3D)*(abs(xd1) + 
     &                     abs(xd2) + abs(xd3))*MACHEPS())
      
            accthr1 = accthr*(dabs(X3(1,1)) + dabs(X3(1,2)) + dabs(X3(1,3)) +
     &                        dabs(X3(2,1)) + dabs(X3(2,2)) + dabs(X3(2,3)) +
     &                        dabs(X3(3,1)) + dabs(X3(3,2)) + dabs(X3(3,3)))/(3d0*dabs(detX3))
     
            if (dabs(real(xd11)).lt.accthr1) xd11 = 0d0
            if (dabs(real(xd12)).lt.accthr1) xd12 = 0d0
            if (dabs(real(xd13)).lt.accthr1) xd13 = 0d0
     
            if (dabs(real(xd22)).lt.accthr1) xd22 = 0d0
            if (dabs(real(xd23)).lt.accthr1) xd23 = 0d0
            if (dabs(real(xd33)).lt.accthr1) xd33 = 0d0
       endif
	   
! End if for alternative D tensor determination	
  	endif

*********************************
* D^{mu nu rho} coefficients 
*	
*********************************

      R01D = 0.5d0*(dc00_1 - dc00_0 - f1D*xd00)
      R02D = 0.5d0*(dc00_2 - dc00_0 - f2D*xd00)
      R03D = 0.5d0*(dc00_3 - dc00_0 - f3D*xd00)

      xd001 = 1d0/detX3*(X3(1,1)*R01D + X3(1,2)*R02D + X3(1,3)*R03D) 	
      xd002 = 1d0/detX3*(X3(2,1)*R01D + X3(2,2)*R02D + X3(2,3)*R03D) 	
      xd003 = 1d0/detX3*(X3(3,1)*R01D + X3(3,2)*R02D + X3(3,3)*R03D) 	

      R11D = 0.5d0*(-dc0_0 - 2d0*dc1_0 - 2d0*dc2_0 - dc11_0 - 2d0*dc12_0 - dc22_0 - f1D*xd11) - 2d0*xd001
      R12D = 0.5d0*(-dc0_0 - 2d0*dc1_0 - 2d0*dc2_0 - dc11_0 - 2d0*dc12_0 - dc22_0 - f2D*xd11) + 0.5d0*dc11_2     
      R13D = 0.5d0*(-dc0_0 - 2d0*dc1_0 - 2d0*dc2_0 - dc11_0 - 2d0*dc12_0 - dc22_0 - f3D*xd11) + 0.5d0*dc11_3
      
      xd111 = 1d0/detX3*(X3(1,1)*R11D + X3(1,2)*R12D + X3(1,3)*R13D)
      xd112 = 1d0/detX3*(X3(2,1)*R11D + X3(2,2)*R12D + X3(2,3)*R13D)
      xd113 = 1d0/detX3*(X3(3,1)*R11D + X3(3,2)*R12D + X3(3,3)*R13D)

      R21D = 0.5d0*(dc1_0 + dc11_0 + dc12_0 - f1D*xd12) - xd002
      R22D = 0.5d0*(dc1_0 + dc11_0 + dc12_0 - f2D*xd12) - xd001
      R23D = 0.5d0*(dc1_0 + dc11_0 + dc12_0 - f3D*xd12) + 0.5d0*dc12_3

c alternative xd112
c      xd112 = 1d0/detX3*(X3(1,1)*R21D + X3(1,2)*R22D + X3(1,3)*R23D)
      xd122 = 1d0/detX3*(X3(2,1)*R21D + X3(2,2)*R22D + X3(2,3)*R23D)
      xd123 = 1d0/detX3*(X3(3,1)*R21D + X3(3,2)*R22D + X3(3,3)*R23D)
      
      R31D = 0.5d0*(dc2_0 + dc12_0 + dc22_0 - f1D*xd13) - xd003
      R32D = 0.5d0*(dc2_0 + dc12_0 + dc22_0 - f2D*xd13) + 0.5d0*dc12_2
      R33D = 0.5d0*(dc2_0 + dc12_0 + dc22_0 - f3D*xd13) - xd001

c alternative xd113, xd123
c      xd113 = 1d0/detX3*(X3(1,1)*R31D + X3(1,2)*R32D + X3(1,3)*R33D)
c      xd123 = 1d0/detX3*(X3(2,1)*R31D + X3(2,2)*R32D + X3(2,3)*R33D)
      xd133 = 1d0/detX3*(X3(3,1)*R31D + X3(3,2)*R32D + X3(3,3)*R33D)
      
      R51D = 0.5d0*(-dc11_0 - f1D*xd22) + 0.5d0*dc11_1
      R52D = 0.5d0*(-dc11_0 - f2D*xd22) - 2d0*xd002
      R53D = 0.5d0*(-dc11_0 - f3D*xd22) + 0.5d0*dc22_3

c alternative xd122
c      xd122 = 1d0/detX3*(X3(1,1)*R51D + X3(1,2)*R52D + X3(1,3)*R53D)
      xd222 = 1d0/detX3*(X3(2,1)*R51D + X3(2,2)*R52D + X3(2,3)*R53D)
      xd223 = 1d0/detX3*(X3(3,1)*R51D + X3(3,2)*R52D + X3(3,3)*R53D)

      R61D = 0.5d0*(-dc12_0 - f1D*xd23) + 0.5d0*dc12_1
      R62D = 0.5d0*(-dc12_0 - f2D*xd23) - xd003
      R63D = 0.5d0*(-dc12_0 - f3D*xd23) - xd002
      
c alternative xd123, xd223
c      xd123 = 1d0/detX3*(X3(1,1)*R61D + X3(1,2)*R62D + X3(1,3)*R63D)
c      xd223 = 1d0/detX3*(X3(2,1)*R61D + X3(2,2)*R62D + X3(2,3)*R63D)
      xd233 = 1d0/detX3*(X3(3,1)*R61D + X3(3,2)*R62D + X3(3,3)*R63D)

      R91D = 0.5d0*(-dc22_0 - f1D*xd33) + 0.5d0*dc22_1
      R92D = 0.5d0*(-dc22_0 - f2D*xd33) + 0.5d0*dc22_2
      R93D = 0.5d0*(-dc22_0 - f3D*xd33) - 2d0*xd003
      
c alternative xd133, xd233
c      xd133 = 1d0/detX3*(X3(1,1)*R91D + X3(1,2)*R92D + X3(1,3)*R93D)
c      xd233 = 1d0/detX3*(X3(2,1)*R91D + X3(2,2)*R92D + X3(2,3)*R93D)
      xd333 = 1d0/detX3*(X3(3,1)*R91D + X3(3,2)*R92D + X3(3,3)*R93D)
      
	end

************************************************************************
      function	A01(q1)
************************************************************************

	double precision q1,q
	double precision A01

#include "../DMNLO_LoopIntegrals.h"

* setting the renorm. scale
	q = muSc

c UV single pole
      if (EpsPole.eq.1d0) then
!--------------------------------------------------------------------            
        A01 = q1**2*UVdiv
c regular contribution
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------
c a0(0)
        if(q1.eq.0.d0)then
          A01 = 0.d0
c a0(m)
        else
          A01 = q1**2*UVdiv + q1**2*(1d0 - dlog(q1**2/q**2))
        endif
!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else
        A01 = 0d0
      endif

      return
      end


************************************************************************
        function  B02(k2,q1,q2)
************************************************************************

	double precision k2,q1,q2,q
	double precision qq2
	complex*16 xp,xm,r1,r2,ieps
	complex*16 B02

#include "../DMNLO_LoopIntegrals.h"

* setting the renorm. scale
	q = muSc
      ieps = dcmplx(0d0,1d-17)

c UV single pole
      if (EpsPole.eq.1d0) then
!--------------------------------------------------------------------            
        B02 = 1d0*UVdiv
c IR single pole
      else if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------            
        if(    (k2.eq.0.d0).and.(q1.eq.0.d0).and.(q2.eq.0.d0)) then
           B02 = -1d0
        else 
           B02 = 0d0
        endif
c regular contribution
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------            
c b0(0;0,0)
      if(    (k2.eq.0.d0).and.(q1.eq.0.d0).and.(q2.eq.0.d0))then
        B02   = UVdiv - IRdiv
c b0(0;0,m) & b0(0;m,0)
      elseif(((k2.eq.0.d0).and.((q1.eq.0.d0).and.(q2.ne.0.d0))).or.
     &       ((k2.eq.0.d0).and.((q2.eq.0.d0).and.(q1.ne.0.d0)))) then
	  qq2 = q1**2 + q2**2
        B02   = UVdiv + 1 - dlog(qq2/q**2)
c b0(0;m,m)
      elseif((k2.eq.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0).and.
     &       (q1.eq.q2  ))then
        B02   = UVdiv - dlog(q1**2/q**2)
c b0(0;q1,q2)
      elseif((k2.eq.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0).and.
     &       (q1.ne.q2  )) then
        B02   = UVdiv + 1d0 - q1**2/(q1**2-q2**2)*dlog(q1**2/q**2)
     &                       -q2**2/(q2**2-q1**2)*dlog(q2**2/q**2)
c b0(k^2;0,0)
      elseif((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.eq.0.d0))then
        B02   = UVdiv + 2d0 - cdlog((-k2-ieps)/q**2)
c b0(m^2;m,0)
      elseif((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.eq.0.d0).and.
     &       (dabs((k2-q1**2)/k2).lt.1.d-4))then
        B02   = UVdiv + 2d0 - dlog(k2/q**2) 
c b0(m^2;0,m)
      elseif((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.ne.0.d0).and.
     &       (dabs((k2-q2**2)/k2).lt.1.d-4))then
        B02   = UVdiv + 2d0 - dlog(k2/q**2)
c b0(k^2;q1,0)
      elseif((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.eq.0.d0).and.
     &       (dabs((k2-q1**2)/k2).gt.1.d-4))then
       B02   = UVdiv + 2d0 - dlog(q1**2/q**2)
     &        - (k2-q1**2)/k2*cdlog((q1**2-k2-ieps)/q1**2)
c b0(k^2;0,q2)
      elseif((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.ne.0.d0).and.
     &       (dabs((k2-q2**2)/k2).gt.1.d-4))then
       B02   = UVdiv + 2d0 - dlog(q2**2/q**2)
     &        - (k2-q2**2)/k2*cdlog((q2**2-k2-ieps)/q2**2)
c b0(k^2;q1,q2)
      elseif((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0)) then

	r1 = 1d0/(2d0*q1*q2)*(-k2+q1**2+q2**2-ieps + 
     &	       cdsqrt(dcmplx((k2-q1**2-q2**2+ieps)**2-4d0*q1**2*q2**2)))
	r2 = 1d0/(2d0*q1*q2)*(-k2+q1**2+q2**2-ieps -
     &	       cdsqrt(dcmplx((k2-q1**2-q2**2+ieps)**2-4d0*q1**2*q2**2)))

	B02 = UVdiv + 2d0 - cdlog(dcmplx(q1*q2/q**2)) + (q1**2-q2**2)/k2*cdlog(dcmplx(q2/q1))
     &           - q1*q2/k2*(r2-r1)*cdlog(r1)		

c Error Message
      else
        write(6,*)'B02 (',sngl(dsqrt(dabs(k2))),',',sngl(q1),',',
     -             sngl(q2),') NOT YET IMPLEMENTED !'
      endif
!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        B02 = 0d0
      endif

      return
      end

************************************************************************
        function  B02p(k2,q1,q2)
************************************************************************

	double precision k2,q1,q2,q
	double precision thmass
	complex*16 xp,xm,kappa,ieps
	complex*16 B02p

#include "../DMNLO_LoopIntegrals.h"

* setting the renorm. scale
	q = muSc
      ieps = dcmplx(0d0,1d-17)

c IR single pole
      if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------            
c b0'(m^2;m,0)
      if((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.eq.0.d0).and.
     &       (dabs((k2-q1**2)/k2).lt.1.d-4)) then
        B02p = -1/(2d0*k2)
        
c b0'(m^2;0,m)
      else if((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.ne.0.d0).and.
     &       (dabs((k2-q2**2)/k2).lt.1.d-4)) then
        B02p = -1/(2d0*k2)
      else
        B02p = 0d0
      endif
c regular contribution
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------            
c b0'(0;0,0)
      if(    (k2.eq.0.d0).and.(q1.eq.0.d0).and.(q2.eq.0.d0)) then
        B02p	= 0.d0
c b0'(0;m,m)
      else if((k2.eq.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0).and.
     &       (q1.eq.q2  )) then
        B02p = 1d0/(6d0*q1**2)
c b0'(0;q1,q2)
      else if((k2.eq.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0).and.
     &       (q1.ne.q2  )) then
        B02p = ((q1**2+q2**2)/2d0 - q1**2*q2**2/(q2**2-q1**2)*
     &		  dlog(q2**2/q1**2))/(q2**2-q1**2)**2
c b0'(k^2;0,0)
      else if((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.eq.0.d0)) then
        B02p = -1d0/k2
c b0'(m^2;m,0)
      elseif((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.eq.0.d0).and.
     &       (dabs((k2-q1**2)/k2).lt.1.d-4)) then
        B02p = 1/(2d0*k2)*(-IRdiv - dlog(q**2)) + 
     &              1/(2d0*k2)*(-2d0 + dlog(k2))

c b0'(m^2;0,m)
      else if((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.ne.0.d0).and.
     &       (dabs((k2-q2**2)/k2).lt.1.d-4)) then
        B02p = 1/(2d0*k2)*(-IRdiv - dlog(q**2)) + 
     &              1/(2d0*k2)*(-2d0 + dlog(k2))

c b0'(k^2;q1,0)
      else if((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.eq.0.d0).and.
     &       (dabs((k2-q1**2)/k2).gt.1.d-4))then
        B02p = 1d0/k2*(-1d0 - q1**2/k2*cdlog((q1**2 - k2 - ieps)/q1**2))

c b0'(k^2;0,q2)
      elseif((k2.ne.0.d0).and.(q1.eq.0.d0).and.(q2.ne.0.d0).and.
     &       (dabs((k2-q2**2)/k2).gt.1.d-4))then
        B02p = 1d0/k2*(-1d0 - q2**2/k2*cdlog((q2**2 - k2 - ieps)/q2**2))

c b0'(k^2;q1,q2)
      elseif((k2.ne.0.d0).and.(q1.ne.0.d0).and.(q2.ne.0.d0))then
        xp  = 1d0/(2d0*k2)*(k2-q2**2+q1**2
     &          + cdsqrt(dcmplx((k2-q2**2+q1**2)**2-4d0*k2*q1**2 + 4d0*k2*ieps)))
        xm  = 1d0/(2d0*k2)*(k2-q2**2+q1**2 
     &          - cdsqrt(dcmplx((k2-q2**2+q1**2)**2-4d0*k2*q1**2 + 4d0*k2*ieps)))

        B02p = 1d0/(k2*(xp-xm))*(xm - xp + xp*(1d0-xp)*cdlog(1d0-1d0/xp) - xm*(1d0-xm)*cdlog(1d0-1d0/xm))

c Error Message
      else
        write(6,*)'B02p(',sngl(dsqrt(dabs(k2))),',',sngl(q1),',',
     -             sngl(q2),') NOT YET IMPLEMENTED !'
      endif
!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        B02p = 0d0
      endif

      return
      end

************************************************************************
      function  C03f(p1,p0,p2,mm0,mm1,mm2)
************************************************************************

      double precision p1,p0,p2,mm0,mm1,mm2		
      double precision k12,s,ss,k22,q1,q2,q3,q,q123		
      double precision old1,old2,old3,oldM1,oldM2,oldM3		

      complex*16 x0,x1,x2,x3,f1,f2,C03reg,C03f,cspen,ckappa, C03
      complex*16 KappaFun,ieps,etad3
      complex*16 t1,t2,FIR1,FIR2,f12,f21,lambda
	  complex*16 zzPat, rzz, szz, ii
	  
      
#include "../DMNLO_LoopIntegrals.h"
		
* setting the renorm. scale
	q = muSc
	ieps = dcmplx(0d0,1d-17)
	ii = dcmplx(0d0,1d0)

	k12 = p1
	s   = p0
	k22 = p2
	q1  = mm0
	q2  = mm1
	q3  = mm2

***************************
* 1. only real part of C0
* 2. IRdiv = 1/eps = 2/(4-D)
* 3. xiIR = 0 - mass regularization of IR (only works for soft poles)

**** argument rotation so that zero masses start at q1

	q123 = q1*q2*q3

	if (q123.eq.0d0) then

	  old1 = k12
	  old2 = s
	  old3 = k22

	  oldM1 = q1
	  oldM2 = q2
	  oldM3 = q3

	 if ((q1.eq.0d0).and.(q2.ne.0d0).and.(q3.eq.0d0)) then

	  k12 = old3
	  k22 = old1
	  q2  = oldM3
	  q3  = oldM2

       else if( (q1.ne.0d0).and.(q2.eq.0d0).and.(q3.eq.0d0) ) then
       
        k12 = old2
        s = old3
        k22 = old1
        q1  = oldM2
        q2  = oldM3	  
        q3  = oldM1	  

	 else if( (q1.ne.0d0).and.(q2.ne.0d0).and.(q3.eq.0d0) ) then

	  k12 = old2
	  s = old1
	  q1  = oldM3
	  q3  = oldM1	  

	 else if( (q1.ne.0d0).and.(q2.eq.0d0).and.(q3.ne.0d0) ) then

	  k22 = old2
	  s = old3
	  q1  = oldM2
	  q2  = oldM1	  

	 endif
	endif
!--------------------------------------------------------------------    
c IR double pole
      if (EpsPole.eq.-2d0) then
!--------------------------------------------------------------------            
c c0(0,s,0;0,0,0)
       if(    (k12.eq.0.d0).and.(k22.eq.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/s
		
c c0(k12,0,0;0,0,0)
       else if(    (k12.ne.0.d0).and.(k22.eq.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k12		
     
c c0(0,0,k22;0,0,0)
       else if(    (k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k22
   
c added from 0308246 soft & collinear 14.5.2013
c c0(0,s,q3^2;0,0,q3)
       else if((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((k22-q3**2)/k22).lt.1.d-4))then
     
        C03f = 1d0/(s-q3**2)*1d0/2d0

c c0(0,q3^2,k22;0,0,q3)
       else if((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((s-q3**2)/s).lt.1.d-4))then
     
        C03f = 1d0/(k22-q3**2)*1d0/2d0
		
       else 
        C03f = 0d0
       endif
!--------------------------------------------------------------------    
c IR single pole
      else if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------    
c c0(0,s,0;0,0,0)
       if(    (k12.eq.0.d0).and.(k22.eq.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/s*cdlog(dcmplx(muSc**2/(-s - ieps)))
     
c c0(k12,0,0;0,0,0)
       else if(    (k12.ne.0.d0).and.(k22.eq.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k12*cdlog(dcmplx(muSc**2/(-k12 - ieps)))
     
c c0(0,0,k22;0,0,0)
       else if(    (k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k22*cdlog(dcmplx(muSc**2/(-k22 - ieps)))
   
c added from 0308246 soft & collinear 14.5.2013
c c0(0,s,q3^2;0,0,q3)
       else if((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((k22-q3**2)/k22).lt.1.d-4))then
     
        C03f = (1d0/(s-q3**2))*(cdlog(muSc**2/(q3**2 - s - ieps)) - (cdlog(dcmplx(muSc**2/(q3**2))))/2d0)

c c0(0,q3^2,k22;0,0,q3)
       else if((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((s-q3**2)/s).lt.1.d-4))then
     
        C03f = (1d0/(k22-q3**2))*(cdlog(muSc**2/(q3**2 - k22 - ieps)) - (cdlog(dcmplx(muSc**2/(q3**2))))/2d0)

c c0(q2^2,s,q3^2;0,q2,q3)
       else if((dabs((k12-q2**2)/k12).lt.1.d-4).and.(dabs((k22-q3**2)/k22).lt.1.d-4).and.(s.ne.0d0).and.
     &            (q1 .eq.0d0).and.(q2.ne.0.d0).and.(q3.ne.0.d0)  )then

        t1 = 1d0/(2d0*s)*(s + q2**2 - q3**2 + cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s)))
        t2 = 1d0/(2d0*s)*(s + q2**2 - q3**2 - cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s)))
        
        lambda = cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s))
        
        FIR1 = 1d0/lambda*(cdlog((s - q2**2 - q3**2 + cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s)))**2/(4d0*q2**2*q3**2)) 
     &        + etad3((t1-1d0)/t1,t2/(t2-1d0),((q2**2-ieps)/s -t2)/((q2**2-ieps)/s -t1)))
        
        C03f = -0.5d0*FIR1
		
c c0(q2^2,0,q2^2;0,q2,q2)
! This special case is needed for DD. Implemented by PS 28.01.2016
       else if((dabs((k12-q2**2)/k12).lt.1.d-4).and.(dabs((k22-q3**2)/k22).lt.1.d-4).and.(s.eq.0d0).and.
     &            (q1.eq.0d0).and.(q2.ne.0.d0).and.(q2.eq.q3))then
        
        C03f = 0.5d0/(q3**2-ieps)
		
c c0(0,s,s;0,0,q3)
! This special case is needed for DD. Implemented by PS 29.01.2016
       else if((k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.
     &       (k22.eq.s).and.(dabs((k22-q3**2)/k22).gt.1.d-4))then
	 	
		zzPat = s/(q3**2 - ieps)
		rzz = 1d0/(1d0 - zzPat)
		C03f = 1d0/q3**2*rzz
		
!		write(*,*)'C03 neu = ',C03f
C 		C03f = C03f + 3d-15
C 		write(*,*)'C03 neu2 = ',C03f
		
!		write(*,*)'Added 0.1 to critical Mandelstam variable to reobtain old case'
C 		s = s + 0.005d0
C 		C03f = 1d0/(s-k22)*cdlog((q3**2 - k22 - ieps)/(q3**2 - s - ieps))
C 		write(*,*)'C0 shift = ', C03f

c implemented from Dittmaier paper    
c c0(0,s,k2^2;0,0,q3)
       else if((k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2 .eq.0.d0).and.(q3 .ne.0.d0).and.
     &       (k22.ne.s).and.(dabs((k22-q3**2)/k22).gt.1.d-4))then
 
        C03f	= 1d0/(s-k22)*cdlog((q3**2 - k22 - ieps)/(q3**2 - s - ieps))
		
       else
        C03f = 0d0
       endif
!--------------------------------------------------------------------    
c regular part
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------               	
	
! DOUBLE POLE functions
c c0(0,s,0;0,0,0)
       if(    (k12.eq.0.d0).and.(k22.eq.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
	 
cc		write(*,*)'Called Special Case C0(0,s,0;0,0,0) with arguments',k12,s,k22,q1,q2,q3
     
        C03f = 1d0/s*(IRdiv**2 + IRdiv*cdlog(dcmplx(muSc**2/(-s - ieps))) 
     &               + (1d0/2d0)*(cdlog(dcmplx(muSc**2/(-s - ieps))))**2 - Lpi**2/6d0) 
		
		
cc		cdlog(dcmplx(muSc**2/(-s - ieps))))**2 - Lpi**2/6d0
		
c c0(k12,0,0;0,0,0)
       else if(    (k12.ne.0.d0).and.(k22.eq.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k12*(IRdiv**2 + IRdiv*cdlog(dcmplx(muSc**2/(-k12 - ieps))) 
     &               + (1d0/2d0)*(cdlog(dcmplx(muSc**2/(-k12 - ieps))))**2 - Lpi**2/6d0) 
	      
c c0(0,0,k22;0,0,0)
       else if(    (k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.eq.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.eq.0.d0)) then
     
        C03f = 1d0/k22*(IRdiv**2 + IRdiv*cdlog(dcmplx(muSc**2/(-k22 - ieps)))
     &              + (1d0/2d0)*(cdlog(dcmplx(muSc**2/(-k22 - ieps))))**2 - Lpi**2/6d0) 
	    
c added from 0308246 soft & collinear 14.5.2013
c c0(0,s,q3^2;0,0,q3)
      elseif((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((k22-q3**2)/k22).lt.1.d-4))then
     
        C03f = (1/(s-q3**2))*((IRdiv**2)/2d0 + IRdiv*(cdlog(dcmplx(muSc**2/(q3**2 - s - ieps))) 
     &        - (cdlog(dcmplx(muSc**2/(q3**2))))/2d0) + ((cdlog(dcmplx(muSc**2/(q3**2 - s - ieps))))**2)/2d0 
     &        - ((cdlog(dcmplx(muSc**2/(q3**2))))**2)/4d0 - cspen(dcmplx((s + ieps)/(s + ieps - q3**2))))

c c0(0,q3^2,k22;0,0,q3)
      elseif((dabs(k12).lt.1.d-4).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.(dabs((s-q3**2)/s).lt.1.d-4))then
     
        C03f = (1/(k22-q3**2))*((IRdiv**2)/2d0 + IRdiv*(cdlog(dcmplx(muSc**2/(q3**2 - k22 - ieps))) 
     &        - (cdlog(dcmplx(muSc**2/(q3**2))))/2d0) + ((cdlog(dcmplx(muSc**2/(q3**2 - k22 - ieps))))**2)/2d0 
     &        - ((cdlog(dcmplx(muSc**2/(q3**2))))**2)/4d0 - cspen(dcmplx((k22 + ieps)/(k22 + ieps - q3**2))))


! SINGLE POLE functions     
! -- soft --

c c0(q2^2,s,q3^2;0,q2,q3)
      else if((dabs((k12-q2**2)/k12).lt.1.d-4).and.(dabs((k22-q3**2)/k22).lt.1.d-4).and.(s.ne.0d0).and.
     &            (q1 .eq.0d0).and.(q2.ne.0.d0).and.(q3.ne.0.d0)  )then
		  
      t1 = 1d0/(2d0*s)*(s + q2**2 - q3**2 + cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s)))
      t2 = 1d0/(2d0*s)*(s + q2**2 - q3**2 - cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s)))
     
      lambda = cdsqrt(dcmplx((s- q2**2 - q3**2)**2 - 4d0*q2**2*q3**2 + 4d0*ieps*s))
     
      FIR1 = 1d0/lambda*(cdlog((s - q2**2 - q3**2 + cdsqrt(dcmplx((s- q2**2 - q3**2)**2 
     &                          - 4d0*q2**2*q3**2 + 4d0*ieps*s)))**2/(4d0*q2**2*q3**2)) 
     &       + etad3((t1-1d0)/t1,t2/(t2-1d0),((q2**2-ieps)/s -t2)/((q2**2-ieps)/s -t1)))
     
      f12 = 0.5d0*cdlog((t1-1d0)/t1)*(cdlog(1d0-t1) + cdlog(-t1) + 2d0*cdlog(t1-t2)) 
     &          - cspen((1d0-t1)/(t2-t1)) + cspen(-t1/(t2-t1))
      f21 = 0.5d0*cdlog((t2-1d0)/t2)*(cdlog(1d0-t2) + cdlog(-t2) + 2d0*cdlog(t2-t1))
     &          - cspen((1d0-t2)/(t1-t2)) + cspen(-t2/(t1-t2))

      FIR2 = FIR1*cdlog(s-ieps) - 1d0/lambda*(f12-f21) 
      
      C03f = -0.5d0*FIR1*IRdiv - 0.5d0*FIR1*cdlog(dcmplx(muSc**2)) + 0.5d0*FIR2
	  
c c0(q2^2,0,q2^2;0,q2,q2)
! This special case is needed for DD. Implemented by PS 28.01.2016
       else if((dabs((k12-q2**2)/k12).lt.1.d-4).and.(dabs((k22-q3**2)/k22).lt.1.d-4).and.(s.eq.0d0).and.
     &            (q1.eq.0d0).and.(q2.ne.0.d0).and.(q2.eq.q3))then
        
	  C03f = 0.5d0*IRdiv/(q3**2-ieps) - 0.5d0/q3**2*cdlog(dcmplx((q3**2-ieps)/muSc**2))

c c0(0,s,s;0,0,q3)
! This special case is needed for DD. Implemented by PS 29.01.2016
       else if((k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1.eq.0.d0).and.(q2.eq.0.d0).and.(q3.ne.0.d0).and.
     &       (k22.eq.s).and.(dabs((k22-q3**2)/k22).gt.1.d-4))then
		
		zzPat = s/(q3**2 - ieps)
 		rzz = 1d0/(1d0 - zzPat)
 		szz = (-ii*Lpi*zzPat + zzPat + zzPat*cdlog(zzPat - 1d0) + cdlog(zzPat - 1d0) - ii*Lpi)/((zzPat - 1d0)*zzPat)
 		C03f = 1d0/q3**2*(IRdiv*rzz + (1d0 + cdlog(muSc**2/(q3**2 - ieps)))*rzz + szz)
!		write(*,*)'C0 neu = ', C03f
		
!		write(*,*)'Added 0.1 to critical Mandelstam variable to reobtain old case'
C 		s = s + 0.1d0
C 		C03f = 1d0/(s-k22)*((IRdiv + cdlog(dcmplx(muSc**2/q3**2)))*cdlog(dcmplx((q3**2 - k22 - ieps)/(q3**2 - s - ieps))) +
C      &       (cdlog(dcmplx((q3**2 - s - ieps)/q3**2)))**2 - (cdlog(dcmplx((q3**2 - k22 - ieps)/q3**2)))**2 +
C      &        cspen(dcmplx((s + ieps)/q3**2)) -  cspen(dcmplx((k22 + ieps)/q3**2)))
C 	 	write(*,*)'C0 shift = ', C03f
		
! -- collinear --
    
c implemented from Dittmaier paper    
c c0(0,s,k2^2;0,0,q3)
      else if((k12.eq.0.d0).and.(k22.ne.0.d0).and.(s.ne.0.d0).and.
     &       (q1 .eq.0.d0).and.(q2 .eq.0.d0).and.(q3 .ne.0.d0).and.
     &       (k22.ne.s).and.(dabs((k22-q3**2)/k22).gt.1.d-4))then
 
        C03f = 1d0/(s-k22)*((IRdiv + cdlog(dcmplx(muSc**2/q3**2)))*cdlog(dcmplx((q3**2 - k22 - ieps)/(q3**2 - s - ieps))) + 
     &       (cdlog(dcmplx((q3**2 - s - ieps)/q3**2)))**2 - (cdlog(dcmplx((q3**2 - k22 - ieps)/q3**2)))**2 + 
     &        cspen(dcmplx((s + ieps)/q3**2)) -  cspen(dcmplx((k22 + ieps)/q3**2)))
		
! CONVERGENT functions (all cases for convergent C0 to be found inside the C03reg function)
      else
           C03f = C03reg(k12,s,k22,q1**2,q2**2,q3**2)
cc		   write (*,*)'Called C03reg with arguments',k12,s,k22,q1**2,q2**2,q3**2
cc		   write (*,*)'Other C0 Call',C03(k12,s,k22,q1,q2,q3)
      endif
!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        C03f = 0d0
      endif
	  
      return
	  
      end

************************************************************************
      function D04f(p1,p21,p32,p3,p2,p13,mm0,mm1,mm2,mm3)
************************************************************************

	  double precision p1,p21,p32,p3,p2,p13,mm0,mm1,mm2,mm3,mm4
	  double precision k1,k21,k32,k3,k2,k13,q0,q1,q2,q3,q4
      double precision k12,k22,k42,s,t,q
      double precision q0123,q123,q23,t1,t2,t3,fac		
	  double precision s_fake, t_fake, u_fake, pcm_fake, costh_fake, kappa

      complex*16 x0,x1,x2,x3,f1,f2,D04,cspen,D04f,ckappa,ieps, LLi2, LLi3
      complex*16 x21,x13,x32,y
      complex*16 KappaFun
	  complex*16 zzPat, rzz, szz

#include "../DMNLO_LoopIntegrals.h"

* setting the renorm. scale
	q = muSc

* small imaginary part
	ieps = dcmplx(0.D0,1.D-17)

cc	write(*,*)'Called D04f with arguments :', p1,p21,p32,p3,p2,p13,mm0,mm1,mm2,mm3

*********************************************
* 1. most of formulas from arXiv:1005.2076
* 2. IRdiv = 2/(4-D) = 1/eps (same as UVdiv)
* 3. we extract a factor (4pi)^eps Gamma(1+eps)

**** argument rotation so that zero masses start at q1
	
	k1  = p1  
	k21 = p21 
	k32 = p32 
	k3  = p3  
	k2  = p2  
	k13 = p13 
	q0  = mm0  
	q1  = mm1  
	q2  = mm2  
	q3  = mm3  	
	
	q0123 = q0*q1*q2*q3

* check if there is a zero mass
	if (q0123.eq.0d0) then
* rotate until first mass is zero
       do while (q0.ne.0d0)
            q4  = q0
            q0  = q1
            q1  = q2  
            q2  = q3  
            q3  = q4  	             
       	
       	t1=k1
       	t2=k2
       	
       	k1  = k21  
       	k21 = k32 
       	k32 = k3 
       	k3  = t1  
       	k2  = k13  
       	k13 = t2 
       enddo

       q123 = q1*q2*q3
* check if there is another zero mass
       if (q123.eq.0d0) then
* rotate until second mass is zero
           do while (q1.ne.0d0)  
             q4  = q1
             q1  = q2
             q2  = q3
             q3  = q4

        	 t1=k1
        	 t2=k21
        	 t3=k3
             
        	 k1  = k2  
        	 k21 = k32 
        	 k32 = k13 
        	 k3  = t1  
        	 k2  = t3  
        	 k13 = t2
        	enddo             
       endif

       q23 = q2*q3
* check if there is another zero mass
       if (q23.eq.0d0) then
* rotate until second mass is zero
           do while (q2.ne.0d0)  
             q4  = q2
             q2  = q3
             q3  = q4  	             

        	 t1=k3
        	 t2=k21

        	 k21 = k13 
        	 k3  = k2  
        	 k2  = t1  
        	 k13 = t2             
        	enddo 
       endif
 	endif

      x13 = -KappaFun(k13+ieps,q1,q3)
      x21 = -KappaFun(k21+ieps,q1,q2)
      x32 = -KappaFun(k32+ieps,q2,q3)

      D04f = 0d0

! 	  write(*,*)'After rotation :', k1,k21,k32,k3,k2,k13,q0,q1,q2,q3
	
!--------------------------------------------------------------------            
c IR DOUBLE POLE
      if (EpsPole.eq.-2d0) then
!--------------------------------------------------------------------               
c soft & collinear, added from Eq (4.25) in 1005.2076, 15.05.2013
c added by JH, QLB
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,m3^2,k2,k13;0,0,m2,m3)
      if ((dabs((k3-q3**2)/k3).lt.1.d-4).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then
     
         D04f = 1d0/((k2 - q2**2)*(k13 - q3**2))
         
c d0(k10,k21,k32,k30,k20,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,m2^2,m3^2;0,0,m2,m3)
c previous case but rotated m0<->m1 (0<->1 indices rotated on k_ij)
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k2-q2**2)/k2).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then

         D04f = (1d0/((k3 - q3**2)*(k21 - q2**2)))

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,m3^2,m3^2,k2,k13;0,0,0,m3)
c added from Eq (4.45) in 1005.2076, 07.08.2021 by LW
      else if ((dabs((k32-q3**2)/k32).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.eq.0d0).and.(k21.eq.0d0).and.
     &        (k2.ne.0d0).and.(k13.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then

         D04f = 2d0/(k2*(k13-q3**2))

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,0,k32,m3^2,0,m3^2;0,0,0,m3)
c added as permutation of zero masses from Eq (4.45) in 1005.2076, 04.03.2022 by LW          
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.ne.0d0).and.(k21.eq.0d0).and.
     &        (k2.eq.0d0).and.(k32.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
       
         D04f = 2d0/(k1*(k32-q3**2))
	 
c soft & collinear, added from Eq (4.29) in 1005.2076, 22.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,k32,k3,k2,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k21.eq.0d0).and.(k2.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 
      D04f = 1d0/(k2*(k13-q3**2))

c added by PS 28.11.14
c previous case but rotated (0<->1)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,0,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k2.eq.0d0).and.(k21.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 
	  D04f = 1d0/(k21*(k3-q3**2))
	 
c soft & collinear, added from Eq (4.19) in 1005.2076, 23.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,m3^2,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
	 		
      D04f = 0.5d0/((k2 -q2**2)*(k13 - q3**2))

c added by PS	5.3.2014
c previous case but rotated (0<->1 and 2<->3)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,k3,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
	 		
      D04f = 0.5d0/((k13 -q3**2)*(k2 - q2**2))

c MISSING SUB-CASES when rotating zero masses among each other...     
     
      else
       D04f = 0d0
      endif

!--------------------------------------------------------------------            
c IR SINGLE POLE
      else if (EpsPole.eq.-1d0) then
!--------------------------------------------------------------------       
c soft, added from Eq (4.4),(4.6) arXiv:1005.2076
c added by KK
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(m1^2,k21,k32,m3^2,k2,k13;0,m1,m2,m3)
      if ((dabs((k1-q1**2)/k1).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &       (q0 .eq.0d0).and.(q1 .ne.0d0).and.(q3 .ne.0d0)) then
		
        if (dabs(dreal(x13)-1d0).ge.1d-13) then

          D04f = dreal(-x13/(q1*q3*(k2 - q2**2)*(1d0 - x13**2))*cdlog(dcmplx(x13)))
c x13 -> 1
        else  

          D04f = 1d0/(2d0*q1*q3*(k2 - q2**2))
          
        endif
		
c added by PS 4.12.14
c previous case but rotated (1<->2)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,k21,k32,m3^2,m2^2,k13;0,m1,m2,m3)
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &       (q0 .eq.0d0).and.(q1 .ne.0d0).and.(q2 .ne.0d0).and.(q3 .ne.0d0)) then
	
		if (dabs(dreal(x32)-1d0).ge.1d-13) then

          D04f = dreal(-x32/(q2*q3*(k1 - q1**2)*(1d0 - x32**2))*cdlog(dcmplx(x32)))

c x32 -> 1
        else  
		
          D04f = 1d0/(2d0*q2*q3*(k1 - q1**2))

       endif
	   
c added by PS 29.2.16
c previous case but rotated (2<->3)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(m1^2,k21,k32,k3,m2^2,k13;0,m1,m2,m3)
      else if ((dabs((k1-q1**2)/k1).lt.1.d-4).and.(dabs((k2-q2**2)/k2).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.ne.0d0).and.(q2.ne.0d0).and.(q3.ne.0d0)) then
	
     	if (dabs(dreal(x21)-1d0).ge.1d-13) then

       	 D04f = dreal(-x21/(q1*q2*(k3 - q3**2)*(1d0 - x21**2))*cdlog(dcmplx(x21)))
		 
c x13 -> 1
     	else  

       	 D04f = 1d0/(2d0*q1*q2*(k3 - q3**2))
      
     	endif

c added by PS	26.02.2016
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,s,s,0,s;0,0,0,m3) 
c This is box B1 for DD in the case of massless quarks -> D0(p1, 0, 0; 0, m1^2, 0,0)
      else if ((k1.eq.0d0).and.(k21.eq.0d0).and.(k32.ne.0d0).and.(dabs((k32 - k3)/k32).lt.1.d-4).and.(k2.eq.0d0).and.(dabs((k32 - k13)/k32).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.eq.0d0 ).and.(q3.ne.0d0)) then
 	
 	    zzPat = k3/(q3**2 - ieps)
 	    rzz = - 1d0/(zzPat - 1d0)**3
 	    D04f = 1d0/q3**4*rzz
				        
c soft, added from Eq (4.6) arXiv:1005.2076
c added by KK
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,k21,k32,m3^2,m1^2,k13;0,0,m1,m3)
c m2 = 0 but not doubly soft divergent as in (4.9)
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &		(dabs((k21-q2**2)/k21).gt.1.d-4).and.(dabs((k13-q3**2)/k13).gt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0).and.(q3 .ne.0d0)) then

          D04f = dreal(-x32/(q2*q3*k1*(1d0 - x32**2))*cdlog(x32))
          
c soft, added from Eq (4.9) arXiv:1005.2076
c added by KK & Saskia
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,m1^2,k32,m3^2,m1^2,m3^2;0,0,m1,m3)
c m2 = 0 doubly soft divergent
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &		(dabs((k21-q2**2)/k21).lt.1.d-4).and.(dabs((k13-q3**2)/k13).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0).and.(q3 .ne.0d0)) then

          D04f = dreal(-2d0*x32/(q2*q3*k1*(1d0 - x32**2))*cdlog(x32))    


c soft & collinear, added from Eq (4.25) in 1005.2076, 15.05.2013
c added by JH, QLB
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,m3^2,k2,k13;0,0,m2,m3)
      else if ((dabs((k3-q3**2)/k3).lt.1.d-4).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then
     
        D04f = (1d0/((k2 - q2**2)*(k13 - q3**2)))*(cdlog((muSc*q2)/(q2**2 - k2 - ieps))
     &      + cdlog((muSc*q3)/(q3**2 - k13 - ieps)))

c d0(k10,k21,k32,k30,k20,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,m2^2,m3^2;0,0,m2,m3)
c previous case but rotated m0<->m1 (0<->1 indices rotated on k_ij)
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k2-q2**2)/k2).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then

        D04f = (1d0/((k3 - q3**2)*(k21 - q2**2)))*(cdlog((muSc*q3)/(q3**2 - k3 - ieps))
     &         + cdlog((muSc*q2)/(q2**2 - k21 - ieps)))

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,m3^2,m3^2,k2,k13;0,0,0,m3)
c added from Eq (4.45) in 1005.2076, 07.08.2021 by LW
      else if ((dabs((k32-q3**2)/k32).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.eq.0d0).and.(k21.eq.0d0).and.
     &        (k2.ne.0d0).and.(k13.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then

        D04f = (1d0/(k2*(k13-q3**2)))*(cdlog(muSc**2/(-k2-ieps)) +2d0*cdlog(muSc*q3/(q3**2-k13-ieps)))

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,0,k32,m3^2,0,m3^2;0,0,0,m3)
c added as permutation of zero masses from Eq (4.45) in 1005.2076, 04.03.2022 by LW         
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.ne.0d0).and.(k21.eq.0d0).and.
     &        (k2.eq.0d0).and.(k32.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
       
         D04f = (1d0/(k1*(k32-q3**2)))*(cdlog(muSc**2/(-k1-ieps)) +2d0*cdlog(muSc*q3/(q3**2-k32-ieps)))

c soft & collinear, added from Eq (4.29) in 1005.2076, 22.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,k32,k3,k2,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k21.eq.0d0).and.(k2.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 
        D04f = 1d0/(k2*(k13-q3**2))*(cdlog(muSc**2/(-k2-ieps)) + cdlog((q3**2-k32-ieps)/(q3**2-k13-ieps)) 
     &		+ cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps)))
	 
c added by PS 28.11.14
c previous case but rotated (0<->1)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,0,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k2.eq.0d0).and.(k21.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 
        D04f = 1d0/(k21*(k3-q3**2))*(cdlog(muSc**2/(-k21-ieps)) + cdlog((q3**2-k32-ieps)/(q3**2-k3-ieps)) 
     &	+ cdlog((q3**2-k13-ieps)/(q3**2-k3-ieps)))
	 
     	 
c soft & collinear, added from Eq (4.19) in 1005.2076, 23.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,m3^2,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
	 
       D04f = 1d0/((k2 -q2**2)*(k13 - q3**2))*(cdlog(muSc*q3/(q3**2-k13-ieps)) 
     &                   + cdlog((q2**2-k21-ieps)/(q2**2-k2-ieps)))

c added by PS	6.3.2014
c previous case but rotated (0<->1 and 2<->3)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,k3,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
	 		
       D04f = 1d0/((k13 -q3**2)*(k2 - q2**2))*(cdlog(muSc*q2/(q2**2-k2-ieps)) 
     &                   + cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps)))
	  
     
c MISSING SUB-CASES when rotating zero masses among each other...     
     
      else
        D04f = 0d0
      endif
     
!--------------------------------------------------------------------            
c REGULAR PART
      else if (EpsPole.eq.0d0) then
!--------------------------------------------------------------------            

c soft, added from Eq (4.4) arXiv:1005.2076, from Eq (4.48) arXiv:0712.1851
c added by KK
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(m1^2,k21,k32,m3^2,k2,k13;0,m1,m2,m3)
      if ((dabs((k1-q1**2)/k1).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.ne.0d0).and.(q2.ne.0d0).and.(q3.ne.0d0)) then

        if (dabs(dreal(x13)-1d0).ge.1d-13) then

          D04f = x13/(q1*q3*(k2 - q2**2)*(1d0 - x13**2))*(
     &      cdlog(dcmplx(x13))*(-IRdiv - dlog(muSc**2))        
     &    - 2d0*cdlog(dcmplx(x13))*cdlog(dcmplx(q2/(q2**2-k2-ieps)))
     &    + cdlog(dcmplx(x21))**2 + cdlog(dcmplx(x32))**2 
     &    - cspen(dcmplx(1d0-x13**2)) 
     &    + LLi3(x13,x21,x32) + LLi3(x13,x21,1d0/x32)
     &    + LLi3(x13,1d0/x21,x32) + LLi3(x13,1d0/x21,1d0/x32))

c x13 -> 1
        else 
		
! This is not necessary for the IR divergent part
		if(x32.eq.x21) then
!			write(*,*)'Added 0.0001 to x32 to reobtain old case'
			x32 = x32 + 0.0001d0
		endif
				 
          D04f = 1d0/(2d0*q1*q3*(k2 - q2**2))*(
     &    (IRdiv + dlog(muSc**2))        
     &    + 2d0*cdlog(dcmplx(q2/(q2**2-k2-ieps) ))
     &    - (1d0 + x21*x32)/(1d0 - x21*x32)*(cdlog(dcmplx(x21)) + cdlog(dcmplx(x32)))
     &    - (x32 + x21)/(x32 - x21)*(cdlog(dcmplx(x21)) - cdlog(dcmplx(x32)))
     &    - 2d0)
	 
       endif

c added by PS 4.12.14
c previous case but rotated (1<->2)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,k21,k32,m3^2,m2^2,k13;0,m1,m2,m3)
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.ne.0d0).and.(q2.ne.0d0).and.(q3.ne.0d0)) then
				
		if (dabs(dreal(x32)-1d0).ge.1d-13) then
		
          D04f = x32/(q2*q3*(k1 - q1**2)*(1d0 - x32**2))*(
     &      cdlog(dcmplx(x32))*(-IRdiv - dlog(muSc**2))        
     &    - 2d0*cdlog(dcmplx(x32))*cdlog(dcmplx(q1/(q1**2-k1-ieps)))
     &    + cdlog(dcmplx(x21))**2 + cdlog(dcmplx(x13))**2 
     &    - cspen(dcmplx(1d0-x32**2)) 
     &    + LLi3(x32,x21,x13) + LLi3(x32,x21,1d0/x13)
     &    + LLi3(x32,1d0/x21,x13) + LLi3(x32,1d0/x21,1d0/x13))
		
c x32 -> 1
        else
			
! This is not necessary for the IR divergent part			
	 	if(x13.eq.x21) then
!			write(*,*)'Added 0.0001 to x13 to reobtain old case'
			x13 = x13 + 0.0001d0
		endif
			
          	D04f = 1d0/(2d0*q2*q3*(k1 - q1**2))*(
     &    (IRdiv + dlog(muSc**2))        
     &    + 2d0*cdlog(dcmplx(q1/(q1**2-k1-ieps) ))
     &    - (1d0 + x21*x13)/(1d0 - x21*x13)*(cdlog(dcmplx(x21)) + cdlog(dcmplx(x13)))
     &    - (x13 + x21)/(x13 - x21)*(cdlog(dcmplx(x21)) - cdlog(dcmplx(x13)))
     &    - 2d0)
	 	 
       endif
	   
c added by PS 29.2.16
c previous case but rotated (2<->3)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(m1^2,k21,k32,k3,m2^2,k13;0,m1,m2,m3)
      else if ((dabs((k1-q1**2)/k1).lt.1.d-4).and.(dabs((k2-q2**2)/k2).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.ne.0d0).and.(q2.ne.0d0).and.(q3.ne.0d0)) then
		
      	if (dabs(dreal(x21)-1d0).ge.1d-13) then

          D04f = x21/(q1*q2*(k3 - q3**2)*(1d0 - x21**2))*(
     &      cdlog(dcmplx(x21))*(-IRdiv - dlog(muSc**2))        
     &    - 2d0*cdlog(dcmplx(x21))*cdlog(dcmplx(q3/(q3**2-k3-ieps)))
     &    + cdlog(dcmplx(x13))**2 + cdlog(dcmplx(x32))**2 
     &    - cspen(dcmplx(1d0-x21**2)) 
     &    + LLi3(x21,x13,x32) + LLi3(x21,x13,1d0/x32)
     &    + LLi3(x21,1d0/x13,x32) + LLi3(x21,1d0/x13,1d0/x32))

c x21 -> 1
       else  
		 
! This is not necessary for the IR divergent part
   	   if(x13.eq.x32) then
!   	   write(*,*)'Added 0.0001 to x32 to reobtain old case'
   		   x32 = x32 + 0.0001d0
   	   endif
		  
         D04f = 1d0/(2d0*q1*q2*(k3 - q3**2))*(
     &    (IRdiv + dlog(muSc**2))        
     &    + 2d0*cdlog(dcmplx(q3/(q3**2-k3-ieps) ))
     &    - (1d0 + x13*x32)/(1d0 - x13*x32)*(cdlog(dcmplx(x13)) + cdlog(dcmplx(x32)))
     &    - (x32 + x13)/(x32 - x13)*(cdlog(dcmplx(x13)) - cdlog(dcmplx(x32)))
     &    - 2d0)
	 
       endif

	
c soft, added from Eq (4.6) arXiv:1005.2076
c added by KK
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,k21,k32,m3^2,m1^2,k13;0,0,m1,m3)
c m2 = 0 but not doubly soft divergent as in (4.9)
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &		(dabs((k21-q2**2)/k21).gt.1.d-4).and.(dabs((k13-q3**2)/k13).gt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0).and.(q3 .ne.0d0)) then

          y = q2*(q3**2-k13-ieps)/(q3*(q2**2-k21-ieps))
     
          D04f = x32/(q2*q3*k1*(1d0 - x32**2))*(
     &      cdlog(x32)*(-IRdiv -0.5d0*cdlog(x32) - cdlog(dcmplx(muSc**2/(q2*q3))) 
     &      - cdlog((q2**2-k21-ieps)/(-k1-ieps)) - cdlog((q3**2-k13-ieps)/(-k1-ieps)))        
     &      + 0.5d0*cdlog(y)**2 - cspen(dcmplx(1d0-x32**2))
     &      + LLi2(x32,y) + LLi2(x32,1d0/y) )
     
c soft, added from Eq (4.9) arXiv:1005.2076
c added by KK & Saskia
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,m1^2,k32,m3^2,m1^2,m3^2;0,0,m1,m3)
c m2 = 0 doubly soft divergent
      else if ((dabs((k2-q2**2)/k2).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.
     &		(dabs((k21-q2**2)/k21).lt.1.d-4).and.(dabs((k13-q3**2)/k13).lt.1.d-4).and.
     &       (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0).and.(q3 .ne.0d0)) then

          D04f = (-2d0*x32/(q2*q3*k1*(1d0 - x32**2))*cdlog(x32))*(IRdiv + cdlog(muSc**2/(-k1-ieps)))          

     
           
c soft & collinear, added from Eq (4.25) in 1005.2076, 15.05.2013
c added by JH, QLB (new version KK 09.07.2013)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,m3^2,k2,k13;0,0,m2,m3)
      else if ((dabs((k3-q3**2)/k3).lt.1.d-4).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then
     
        D04f = (1d0/((k2 - q2**2)*(k13 - q3**2)))*(IRdiv**2 + IRdiv*(cdlog((muSc*q2)/(q2**2 - k2 - ieps))
     &      + cdlog((muSc*q3)/(q3**2 - k13 - ieps))) 
     &      + 2d0*cdlog((muSc*q2)/(q2**2 - k2 - ieps))*cdlog((muSc*q3)/(q3**2 - k13 - ieps))
     &      - cdlog(x32)**2 - 2d0*Lpi**2/3d0)

c d0(k10,k21,k32,k30,k20,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,m2^2,m3^2;0,0,m2,m3)
c previous case but rotated m0<->m1 (0<->1 indices rotated on k_ij)
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k2-q2**2)/k2).lt.1.d-4).and.(k1.eq.0d0).and.
     &        (q0.eq.0d0).and.(q1.eq.0d0).and.(q2.ne.0d0 ).and.(q3.ne.0d0)) then
		
        D04f = (1d0/((k3 - q3**2)*(k21 - q2**2)))*(IRdiv**2 + IRdiv*(cdlog((muSc*q3)/(q3**2 - k3 - ieps))
     &      + cdlog((muSc*q2)/(q2**2 - k21 - ieps))) 
     &      + 2d0*cdlog((muSc*q3)/(q3**2 - k3 - ieps))*cdlog((muSc*q2)/(q2**2 - k21 - ieps))
     &      - cdlog(x32)**2 - 2d0*Lpi**2/3d0)

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,m3^2,m3^2,k2,k13;0,0,0,m3)
c added from Eq (4.45) in 1005.2076, 07.08.2021 by LW
      else if ((dabs((k32-q3**2)/k32).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.eq.0d0).and.(k21.eq.0d0).and.
     &        (k2.ne.0d0).and.(k13.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then

        D04f = (1d0/(k2*(k13-q3**2)))*(2d0*IRdiv**2+IRdiv*(cdlog(muSc**2/(-k2-ieps)) +2d0*cdlog(muSc*q3/(q3**2-k13-ieps)))
     &      + 2d0*cdlog(muSc**2/(-k2-ieps))*cdlog(muSc*q3/(q3**2-k13-ieps))-5d0*Lpi**2/6d0)

c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(k1,0,k32,m3^2,0,m3^2;0,0,0,m3)
c added as permutation of zero masses from Eq (4.45) in 1005.2076, 04.03.2022 by LW         
      else if ((dabs((k13-q3**2)/k13).lt.1.d-4).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(k1.ne.0d0).and.(k21.eq.0d0).and.
     &        (k2.eq.0d0).and.(k32.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
       
         D04f = (1d0/(k1*(k32-q3**2)))*(2d0*IRdiv**2+IRdiv*(cdlog(muSc**2/(-k1-ieps)) +2d0*cdlog(muSc*q3/(q3**2-k32-ieps)))
     &      + 2d0*cdlog(muSc**2/(-k1-ieps))*cdlog(muSc*q3/(q3**2-k32-ieps))-5d0*Lpi**2/6d0)
	 
c soft & collinear, added from Eq (4.29) in 1005.2076, 22.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,k32,k3,k2,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k21.eq.0d0).and.(k2.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 	 
        D04f = 1d0/(k2*(k13-q3**2))*(IRdiv**2 + IRdiv*(cdlog(muSc**2/(-k2-ieps)) + cdlog((q3**2-k32-ieps)/(q3**2-k13-ieps)) 
     &	+ cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps))) + 0.5d0*(cdlog(muSc**2/(-k2-ieps)))**2
     &	+ cdlog(muSc**2/(-k2-ieps))*( cdlog((q3**2-k32-ieps)/(q3**2-k13-ieps)) + cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps))) 
     &	- 0.5d0*(cdlog((q3**2-k32-ieps)/(q3**2-k3-ieps)))**2 
     &      - 2d0*cspen((k32-k13)/(q3**2-k13-ieps)) - 2d0*cspen((k3-k13)/(q3**2-k13-ieps))
     &	+ LLi2((k2-ieps)/(k32-ieps-q3**2),q3**2/(q3**2-k3-ieps)) - Lpi**2/3d0)
     
! + Lpi**2/6d0)    for comparison with LoopTools

c added by PS 28.11.14
c previous case but rotated (0<->1)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,k3,0,k13;0,0,0,m3)
      else if ((k1.eq.0d0).and.(k2.eq.0d0).and.(k21.ne.0d0).and.(q0.eq.0d0).and.(q1.eq.0d0).and.(q2.eq.0d0 ).and.(q3.ne.0d0)) then
	 	 
        D04f = 1d0/(k21*(k3-q3**2))*(IRdiv**2 + IRdiv*(cdlog(muSc**2/(-k21-ieps)) + cdlog((q3**2-k32-ieps)/(q3**2-k3-ieps)) 
     &	+ cdlog((q3**2-k13-ieps)/(q3**2-k3-ieps))) + 0.5d0*(cdlog(muSc**2/(-k21-ieps)))**2
     &	+ cdlog(muSc**2/(-k21-ieps))*( cdlog((q3**2-k32-ieps)/(q3**2-k3-ieps)) + cdlog((q3**2-k13-ieps)/(q3**2-k3-ieps))) 
     &	- 0.5d0*(cdlog((q3**2-k32-ieps)/(q3**2-k13-ieps)))**2 
     &      - 2d0*cspen((k32-k3)/(q3**2-k3-ieps)) - 2d0*cspen((k13-k3)/(q3**2-k3-ieps))
     &	+ LLi2((k21-ieps)/(k32-ieps-q3**2),q3**2/(q3**2-k13-ieps)) - Lpi**2/3d0)
     
! + Lpi**2/6d0)    for comparison with LoopTools

c added by PS	26.02.2016
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,0,s,s,0,s;0,0,0,m3) 
c This is box B1 for DD in the case of massless quarks -> D0(p1, 0, 0; 0, m1^2, 0,0)
      else if ((k1.eq.0d0).and.(k21.eq.0d0).and.(k32.ne.0d0).and.(dabs((k32 - k3)/k32).lt.1.d-4).and.(k2.eq.0d0).and.(dabs((k32 - k13)/k32).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.eq.0d0 ).and.(q3.ne.0d0)) then

	 	zzPat = k3/(q3**2 - ieps)
	 	rzz = - 1d0/(zzPat-1)**3
	 	szz = (2d0*cdlog(1d0 - zzPat) - zzPat)/(zzPat - 1d0)**3
	 	D04f = 1d0/q3**4*(IRdiv*rzz + (1d0 + cdlog(muSc**2/(q3**2 - ieps)))*rzz + szz)
			 
c soft & collinear, added from Eq (4.19) in 1005.2076, 23.05.2013
c added by PS
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,k21,k32,m3^2,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k3-q3**2)/k3).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
	 		
        D04f = 1d0/((k2 -q2**2)*(k13 - q3**2))*(0.5d0*IRdiv**2 + IRdiv*(cdlog(muSc*q3/(q3**2-k13-ieps)) 
     &      + cdlog((q2**2-k21-ieps)/(q2**2-k2-ieps))) 
     &	+ 2d0*cdlog((q2**2-k21-ieps)/(q2**2-k2-ieps))*cdlog(muSc*q3/(q3**2-k13-ieps)) 
     &      + (cdlog(muSc*q3/(q3**2-k13-ieps)))**2 - 2d0*cspen((k21-k2)/(q2**2-k2-ieps))
     &	+ LLi2((q3**2-k13-ieps)/(q2**2-k21-ieps),(q2/q3)*x32) + LLi2((q3**2-k13-ieps)/(q2**2-k21-ieps),(q2/q3)/x32)
     &	- Lpi**2/6d0)
     
! + Lpi**2/12.0d0) for comparison with LoopTools (1/12 instead of 1/6 due to the factor 1/2 in front of IRdiv**2)
     
c added by PS	6.3.2014
c previous case but rotated (0<->1 and 2<->3)
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3) = d0(0,m2^2,k32,k3,k2,k13;0,0,m2,m3)
      else if ((k1.eq.0d0).and.(dabs((k21-q2**2)/k21).lt.1.d-4).and.(q0.eq.0d0).and.(q1.eq.0d0).and.
     &        (q2.ne.0d0 ).and.(q3.ne.0d0)) then
		
		D04f = 1d0/((k13 -q3**2)*(k2 - q2**2))*(0.5d0*IRdiv**2 + IRdiv*(cdlog(muSc*q2/(q2**2-k2-ieps)) 
     &      + cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps))) 
     &	+ 2d0*cdlog((q3**2-k3-ieps)/(q3**2-k13-ieps))*cdlog(muSc*q2/(q2**2-k2-ieps)) 
     &      + (cdlog(muSc*q2/(q2**2-k2-ieps)))**2 - 2d0*cspen((k3-k13)/(q3**2-k13-ieps))
     &	+ LLi2((q2**2-k2-ieps)/(q3**2-k3-ieps),(q3/q2)*x32) + LLi2((q2**2-k2-ieps)/(q3**2-k3-ieps),(q3/q2)/x32)
     &	- Lpi**2/6d0)
	 
! + Lpi**2/12.0d0) for comparison with LoopTools (1/12 instead of 1/6 due to the factor 1/2 in front of IRdiv**2)
	 
c MISSING SUB-CASES when rotating zero masses among each other...     

c regular case     
c d0(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3)
c      else if ((q0 .ne.0d0).and.(q1 .ne.0d0).and.(q2 .ne.0d0).and.(q3 .ne.0d0)) then

c WARNING this can cause problems for unimplemented IR-divergent cases !!!!!!!!!!!!!!!!!!!!!!

      else 
	  
c		write(*,*)'Called Denner D0 with arguments', k1,k21,k32,k3,k2,k13,q0,q1,q2,q3
c		write(*,*)'Called Denner D0'
						
! These cases occur in B2/B4 in DD, we add a little bit of noise to avoid numerical instabilities in the Denner function 
		if ((k1.eq.q3**2).and.(k21.eq.0d0).and.(k32.eq.k13).and.(k2.eq.q3**2).and.(q0.ne.0d0).and.(q1.ne.0d0).and.(q2.ne.0d0).and.(q3.ne.0d0)) then
			
			pcm_fake = 0.5
			costh_fake = 0.5
			
			s_fake = k13 + 2.0d0*pcm_fake**2 + q3**2 + 2.0d0*dsqrt(k13 + pcm_fake**2)*dsqrt(q3**2 + pcm_fake**2)
			t_fake = 2.0d0*k13 - 0.5d0*(s_fake + k13 - q3**2)**2/s_fake + 0.5d0*kappa(s_fake, k13, q3**2)**2/s_fake*costh_fake
			u_fake = 2.0d0*k13**2 + 2.0d0*q3**2 - t_fake - s_fake
			
			if(abs(k3 - (dsqrt(k13) + q3)**2).lt.0.01d0) then
				write(*,*)'Added some small numerical noise to t & s to improve numerical stability of Denner D0.'
				k3 = s_fake
				k21 = t_fake
			endif
			
			if(abs(k3 - (dsqrt(k13) - q3)**2).lt.0.01d0) then
				write(*,*)'Added some small numerical noise to t & u to improve numerical stability of Denner D0.'
				k3 = u_fake
				k21 = t_fake
			endif	
		endif
	
        D04f = D04(k1,k21,k32,k3,k2,k13,q0,q1,q2,q3)

		
	endif

!--------------------------------------------------------------------
c all other pole coefficients which are zero
      else 
        D04f = 0d0
      endif


	return
	end

************************************************************************
        FUNCTION C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
* Es muss die Dittmaier-Routine D04 hinzugelinkt werden fuer die eta   *
*                   und die Spence-Funktionen                          *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
c ----------------------------------------------
      complex*16 cx,cy
c ----------------------------------------------
      COMPLEX*16 C03,CSPEN,ETAD,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM
C     REAL*8 KAPPA
      COMPLEX*16 KAPPA
C     KAPPA(A,B,C) = DSQRT(A**2+B**2+C**2-2*(A*B+A*C+B*C))
C     KAPPA(A,B,C) = DSQRT(DABS(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      KAPPA(A,B,C) = CDSQRT(DCMPLX(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      EPS = 1.D-8
      IM = DCMPLX(0.D0,1.D0)
      IEPS = DCMPLX(0.D0,1.D-17)
      PI = 4*DATAN(1.D0)
      XX = 0.D0
C     IF(P1.LT.0.D0.OR.P2.LT.0.D0.OR.P3.LT.0.D0) XX=1.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
c ----------------------------------------------
      alpha = cdabs( kappa(q10,q21,q20) )
c      ALPHA = KAPPA(Q10,Q21,Q20)
c ----------------------------------------------
      ALP(0) = KAPPA(Q21,SM1,SM2)*(1+IEPS*Q21)
      ALP(1) = KAPPA(Q20,SM2,SM0)*(1+IEPS*Q20)
      ALP(2) = KAPPA(Q10,SM0,SM1)*(1+IEPS*Q10)
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
c ----------------------------------------------
      cdum = dcmplx(0.D0)
c      CDUM=0.D0
c ----------------------------------------------
      DO I=0,2
       DO J=1,2
        CDUM = CDUM + CSPEN((Y0(I)-1)/Y(I,J)) - CSPEN(Y0(I)/Y(I,J))
        CX = ETAD(1-X(I,J),1/Y(I,J))
c ----------------------------------------------
        if (cdabs(cx).ne.0.D0) then 
c        IF(CX.NE.0.D0)THEN
c ----------------------------------------------
         CDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = ETAD(-X(I,J),1/Y(I,J))
c ----------------------------------------------
        if (cdabs(cy).ne.0.D0) then 
c        IF(CY.NE.0.D0)THEN 
c ----------------------------------------------
         CDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = ETAD(-X(I,1),-X(I,2))
c ----------------------------------------------
       if (cdabs(cx).ne.0.D0) then 
c       IF(CX.NE.0.D0)THEN
c ----------------------------------------------
        CDUM = CDUM - CX*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       CY = ETAD(Y(I,1),Y(I,2))
c ----------------------------------------------
       if (cdabs(cy).ne.0.D0) then 
c       IF(CY.NE.0.D0)THEN
c ----------------------------------------------
        CDUM = CDUM + CY*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2*PI*IM*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      C03 = CDUM/ALPHA
      RETURN
      END
c ======================================================================

************************************************************************
        FUNCTION C03reg(P1,P0,P2,M0,M1,M2)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P0,P1,P2 = SQUARED EXTERNAL MOMENTA  			                 *
*  M0,M1,M2 = SQUARED LOOP MASSES          			           *
*----------------------------------------------------------------------*
* based on the function from 't Hooft Veltman (& Hitchhiker's guide)
* (only for real masses and momenta)
*----------------------------------------------------------------------*
* 29.5.2013   K.Kovarik                                                *
************************************************************************
      Implicit None
      REAL*8 P0,P1,P2,M0,M1,M2
c ----------------------------------------------
      REAL*8 Pt0,Pt1,Pt2,Mt0,Mt1,Mt2
      COMPLEX*16 C03reg,CSPEN,ETAD
      COMPLEX*16 CDUM,cI,kk
      complex*16 aa,aa1,aa2,dd,y0temp,y0,y0Den0,y0Den1,y0Den2,yp,ym
      complex*16 a,b,c,xm,xp
      double precision ak,bk,ck,eps,MACHEPS
      complex*16 temp,y0new
      complex*16 fa,fb,fc,fd,fe,ff
      COMPLEX*16 KAPPA,S3int,S2int
	  
      KAPPA(ak,bk,ck) = CDSQRT(DCMPLX(ak**2+bk**2+ck**2-2*(ak*bk+ak*ck+bk*ck)))

      eps= MACHEPS()
      cI =(0d0,1d0)
	  
cc	  write(*,*)'M0 =',M0
cc	  write(*,*)'M1 =',M1
cc	  write(*,*)'M2 =',M2
		
*********************** checking for zero squared momenta
      if (P0*P1*P2 .eq. 0d0) then
c saving original arguments
            Pt0=P0  
            Pt1=P1  
            Pt2=P2  
            Mt0=M0  
            Mt1=M1  
            Mt2=M2
        if (P0.eq.0d0 .and. P1.eq.0d0 .and. P2.eq.0d0) then              
*********************** C0 with all zero squared momenta
cc			write(*,*)'All three momenta zero'
              if (M0.eq.M1) then
                    if (M1.eq.M2) then
					  if (M2.ne.0d0) then
cc						write(*,*)'CASE 9'
						C03reg = -1d0/M2
					  else
!						write(*,*)'WARNING! Encountered C0 with all arguments zero. This has to be discussed in terms of IR & UV pole cancellation.'
						C03reg = 0d0
					  endif     
                    else
cc					  write(*,*)'CASE 8'
                      C03reg = 1d0/(M2-M1)*(1d0 - M2/(M1-M2)*dlog(M1/M2))
                    endif
              else if (M1.eq.M2) then
cc					 write(*,*)'CASE 8'
                     C03reg = 1d0/(M1-M0)*(-1d0+ M0/(M0-M1)*dlog(M0/M2))
              else 
cc				write(*,*)'CASE 7'
                C03reg = 1d0/(M1-M0)*((M0*dlog(M0)-M2*dlog(M2))/(M0-M2) - (M1*dlog(M1)-M2*dlog(M2))/(M1-M2))
              endif

        else if (P0*P1 + P0*P2 + P2*P1 .eq. 0d0) then
*********************** C0 with 1 non-zero squared momenta (rotating so that P1=0,P0=0)
            if (P0.eq.0d0.and.P2.eq.0d0) then
c rotated P2<->P1 & M1<->M2
              Pt1=P2  
              Pt2=P1  
              Mt1=M2  
              Mt2=M1            
            else if (P1.eq.0d0.and.P2.eq.0d0) then
c rotated P2<->P0 & M1<->M0
              Pt0=P2  
              Pt2=P0  
              Mt0=M1  
              Mt1=M0            
            endif
         
c--- C0 parameters      
c-- tHooft Veltman/Denner habilitation choice (Hitchhiker's guide Tab.2.1)
         fa = dcmplx(Pt0)
         fb = dcmplx(Pt1)
         fc = dcmplx(Pt2-Pt1-Pt0)
         fd = dcmplx(Mt1-Pt0-Mt2)
         fe = dcmplx(Pt0-Pt2-Mt1+Mt0)
         ff = dcmplx(Mt2)-cI*eps**2

c---------------------
         a=fc
         b=fd+fe
         c=ff
         
         y0 = -fe/fc
c only for Denner choice of a,b,c,d,e,f      
         yp = (-Mt0 + Mt2 + Pt2 + CDSQRT(DCMPLX(-4*Mt2*Pt2 + (-Mt0 + Mt2 + Pt2)**2 + cI*eps**2*Pt2)))/(2.*Pt2)
         ym = (-Mt0 + Mt2 + Pt2 - CDSQRT(DCMPLX(-4*Mt2*Pt2 + (-Mt0 + Mt2 + Pt2)**2 + cI*eps**2*Pt2)))/(2.*Pt2)
         
         temp = S3int(a,b,c,y0,yp,ym)
         C03reg = -1d0/fc*temp
         

c cases C0(0,0,p,m,m,m) & C0(0,0,p,m_0,m,m)            
         if (Mt1.ne.Mt2) then 
           temp = S2int(fd,ff,y0)
           C03reg = C03reg + 1d0/fc*temp
         endif
         
            
        else
*********************** C0 with 2 non-zero squared momenta (rotating so that P1=0)
cc	write(*,*)'ONE ZERO SQUARED MOMENTA'
        if (P0.eq.0d0) then
c rotated P0<->P1 & M0<->M2
            Pt0=P1  
            Pt1=P0  
            Mt0=M2  
            Mt2=M0            
        else if (P2.eq.0d0) then
c rotated P2<->P1 & M1<->M2
            Pt2=P1  
            Pt1=P2  
            Mt1=M2  
            Mt2=M1            
        endif
cc		write(*,*)'One momentum zero'
		 
c--- C0 parameters      
c-- tHooft Veltman/Denner habilitation choice (Hitchhiker's guide Tab.2.1)
        fa = dcmplx(Pt0)
        fb = dcmplx(Pt1)
        fc = dcmplx(Pt2-Pt1-Pt0)
        fd = dcmplx(Mt1-Pt0-Mt2)
        fe = dcmplx(Pt0-Pt2-Mt1+Mt0)
        ff = dcmplx(Mt2)-cI*eps**2
		
      if (Pt0.eq.Pt2) then
            if (Mt1.eq.Mt0) then
c---------------------
cc			write(*,*)'CASE 3'
            a=fa
            b=fd
            c=ff
			
			temp = cdsqrt(4*a*c-b**2)
					
			C03reg = -1d0/(2*a)*cdlog(a+b+c) + 1d0/(2*a)*cdlog(c) + b/(a*temp)*atan((2*a+b)/temp) - b/(a*temp)*atan(b/temp)
			
            else      
c---------------------
cc			write(*,*)'CASE 2'
            a=fa
            b=fd+fe
            c=ff

         xp = (-b+cdsqrt(b**2-4*a*c))/(2*a)
         xm = (-b-cdsqrt(b**2-4*a*c))/(2*a)

         temp = cdlog((1d0-xp)*(1d0-xm)) - xp*cdlog(1d0-1d0/xp) - xm*cdlog(1d0-1d0/xm)

         C03reg = - 1d0/fe*temp

c---------------------
            a=fa
            b=fd
            c=ff

         xp = (-b+cdsqrt(b**2-4*a*c))/(2*a)
         xm = (-b-cdsqrt(b**2-4*a*c))/(2*a)

         temp = cdlog((1d0-xp)*(1d0-xm)) - xp*cdlog(1d0-1d0/xp) - xm*cdlog(1d0-1d0/xm)

         C03reg = C03reg + 1d0/fe*temp

      endif
      else

                
c---------------------
         a=fa
         b=fd
         c=ff
         
         y0 = -fe/fc
c only for Denner choice of a,b,c,d,e,f      
         yp = (-Mt1 + Mt2 + Pt0 + CDSQRT(DCMPLX(-4*Mt2*Pt0 + (-Mt1 + Mt2 + Pt0)**2 + cI*eps**2*Pt0)))/(2.*Pt0)
         ym = (-Mt1 + Mt2 + Pt0 - CDSQRT(DCMPLX(-4*Mt2*Pt0 + (-Mt1 + Mt2 + Pt0)**2 + cI*eps**2*Pt0)))/(2.*Pt0)
         
         temp = S3int(a,b,c,y0,yp,ym)
         C03reg = 1d0/fc*temp

c---------------------
         a=fa+fc
         b=fd+fe
         c=ff
         
         y0 = -fe/fc
c only for Denner choice of a,b,c,d,e,f      
         yp = (-Mt0 + Mt2 + Pt2 + CDSQRT(DCMPLX(-4*Mt2*Pt2 + (-Mt0 + Mt2 + Pt2)**2 + cI*eps**2*Pt2)))/(2.*Pt2)
         ym = (-Mt0 + Mt2 + Pt2 - CDSQRT(DCMPLX(-4*Mt2*Pt2 + (-Mt0 + Mt2 + Pt2)**2 + cI*eps**2*Pt2)))/(2.*Pt2)
         
         temp = S3int(a,b,c,y0,yp,ym)
         C03reg = C03reg - 1d0/fc*temp
                    
        endif
      endif
            
      else
*********************** C0 with no zero squared momenta
cc	write(*,*)'NO ZERO SQUARED MOMENTA'
c--- C0 parameters      
c-- tHooft Veltman/Denner habilitation choice (Hitchhiker's guide Tab.2.1)
      fa = dcmplx(P0)
      fb = dcmplx(P1)
      fc = dcmplx(P2-P1-P0)
      fd = dcmplx(M1-P0-M2)
      fe = dcmplx(P0-P2-M1+M0)
      ff = dcmplx(M2)-cI*eps**2
            
c----------      
      kk    = Kappa(P0,P1,P2)
	  
!	  write(*,*)'Kappa = ',kk
	  
! PS: Somewhat dirty, but efficient hack to cure problems for Kappa = 0.
	  if(abs(kk).lt.0.1d0) then
!		write(*,*)'WARNING! |Kappa| < 0.1, automatically set to 0.1. Old Kappa = ', kk
		kk = 0.1d0
	  endif
	  
c----------      
      aa = 1d0/(2d0*P1)*(P0+P1-P2+kk)

! Check if Kappa = 0, otherwise y0 etc. is not defined	  
!	  if(abs(kk).ne.0d0) then
	  
      y0temp = (-fd-fe*aa)/(fc+2*fb*aa)  
      
cc	  write(*,*)'y0 =',y0temp 
c---------------------
      a=fb
      b=fc+fe
      c=fa+fd+ff
      y0=y0temp+aa
      
cc	  write(*,*)'y00 = ',y0
	  
c only for Denner choice of a,b,c,d,e,f      
      y0 = -(M0*P0 + M0*P1 - 2*M2*P1 + P0*P1 - P1**2 - kk*(-M0 + M1 + P1) - M0*P2 + 
     &     P1*P2 + M1*(-P0 + P1 + P2))/(2.*kk*P1)
      yp = (-M0 + M1 + P1 + CDSQRT(DCMPLX(-4*M1*P1 + (-M0 + M1 + P1)**2 + cI*eps**2*P1)))/(2.*P1)
      ym = (-M0 + M1 + P1 - CDSQRT(DCMPLX(-4*M1*P1 + (-M0 + M1 + P1)**2 + cI*eps**2*P1)))/(2.*P1)
      
      
      temp = S3int(a,b,c,y0,yp,ym)      
      C03reg = (-1d0/kk)*temp

c--------------
      a=fa+fb+fc
      b=fd+fe
      c=ff
      y0=y0temp/(1d0-aa)
      
cc	  write(*,*)'y01 = ',y0

c only for Denner choice of a,b,c,d,e,f      
      y0 = (M0*P0 - M2*P0 - M0*P1 + M2*P1 + M0*P2 - 2*M1*P2 + M2*P2 + P0*P2 + P1*P2 - 
     &    P2**2 + kk*(-M0 + M2 + P2))/(2.*kk*P2)
      yp = (-M0 + M2 + P2 + CDSQRT(DCMPLX(-4*M2*P2 + (-M0 + M2 + P2)**2 + cI*eps**2*P2)))/(2.*P2)
      ym = (-M0 + M2 + P2 - CDSQRT(DCMPLX(-4*M2*P2 + (-M0 + M2 + P2)**2 + cI*eps**2*P2)))/(2.*P2)
      
      temp = S3int(a,b,c,y0,yp,ym)
      C03reg = C03reg - (-1d0/kk)*temp

c--------------
      a=fa
      b=fd
      c=ff
      y0=-y0temp/aa
	 
cc	  write(*,*)'y02 = ',y0

c only for Denner choice of a,b,c,d,e,f      
      y0 = (2*M0*P0 - M2*P0 + P0**2 + kk*(-M1 + M2 + P0) - (M2 + P0)*P1 + M2*P2 - 
     &    P0*P2 - M1*(P0 - P1 + P2))/(2.*kk*P0)
      yp = (-M1 + M2 + P0 + CDSQRT(DCMPLX(-4*M2*P0 + (-M1 + M2 + P0)**2 + cI*eps**2*P0)))/(2.*P0)
      ym = (-M1 + M2 + P0 - CDSQRT(DCMPLX(-4*M2*P0 + (-M1 + M2 + P0)**2 + cI*eps**2*P0)))/(2.*P0)

      temp = S3int(a,b,c,y0,yp,ym)
      C03reg = C03reg + (-1d0/kk)*temp

      endif

      END

!***************************************************************************
* Function to compute (according to 't Hooft Veltman or Hitchhiker's guide)
* 
*    int( ( ln(a*y^2+b*y+c) - ln(a*y0^2+b*y0+c) )/(y-y0) ,y=0..1 )
*
!***************************************************************************

      FUNCTION S3int(a,b,c,y0,yp,ym)

      COMPLEX*16 S3int,a,b,c,y0,res,reseta,reseta2,yp,ym
      COMPLEX*16 y1,y2,d,cI,zz
      COMPLEX*16 CSPEN,ETAD,ETAD3
      double precision epsn,MACHEPS,pi
      
      external Qsol

      cI  = (0d0,1d0)
      pi = 4*DATAN(1D0)
      epsn = MACHEPS()

      y1 = yp
      y2 = ym

        res = CSPEN(y0/(y0-y1)) + CSPEN(y0/(y0-y2)) - CSPEN((y0-DCMPLX(1D0))/(y0-y1))
     &        - CSPEN((y0-DCMPLX(1D0))/(y0-y2))
        
        reseta = (0d0,0d0)
        if (dreal(y0).ne.1d0) then
        
        reseta = ETAD(-y1,1d0/(y0-y1))*cdlog(y0/(y0-y1)) 
     &         - ETAD(DCMPLX(1D0)-y1,1d0/(y0-y1))*cdlog((y0-DCMPLX(1D0))/(y0-y1))
     &         + ETAD(-y2,1d0/(y0-y2))*cdlog(y0/(y0-y2)) 
     &         - ETAD(DCMPLX(1D0)-y2,1d0/(y0-y2))*cdlog((y0-DCMPLX(1D0))/(y0-y2))
                
c------ Cancellations in complex parts (can spoil the sign of the small ieps complex part)
c
c a*(y0-y_i^+)*(y0-y_i^-) = a*y_0^2 + b*y_0 + c

        zz = a*y0**2+b*y0
        zzreal = dabs(dble(0.5d0*(zz+dconjg(zz))))
        zzimg  = aimag(0.5d0*(zz-dconjg(zz)))
        
c Cancellations are of the order of epsn - machine epsilon - the smallest positive number 
c    if the imaginary part of zz = a*y0**2+b*y0 is due to cancellations erase it
        if (zzimg.le.zzreal*epsn*1d2 .and. zzimg.ge.zzreal*epsn/1d4) then
              zz = zzreal + c
        else
              zz = zz + c
        endif
                
        reseta2 = ETAD3(-y1,-y2,c/a) - ETAD3(y0-y1,y0-y2,zz)
        if (dreal(a).le.0d0 .and. aimag(zz).ge.0d0) then
              reseta2 = reseta2 - 2*pi*cI
        endif
         
        reseta = reseta + reseta2*cdlog((y0-DCMPLX(1D0))/y0)

        endif

        S3int = res + reseta
                
      END
      
!***************************************************************************
* Function to compute (according to 't Hooft Veltman or Hitchhiker's guide)
*    int( ( ln(a*y+b) - ln(a*y0+b) )/(y-y0) ,y=0..1 )
!***************************************************************************

      FUNCTION S2int(a,b,y0)

      COMPLEX*16 S2int,a,b,y0,res,reseta
      COMPLEX*16 y1,cI,zz
      COMPLEX*16 CSPEN,ETAD,ETAD3
      double precision epsn,MACHEPS,pi

      external Qsol

      cI  = (0d0,1d0)
      pi = 4*DATAN(1D0)
      epsn = MACHEPS()

      if (dreal(a).eq.0d0) then

        S2int = cdlog(b)*cdlog((y0-DCMPLX(1D0))/y0)      
      
      else
        y1 = -b/a

        res = CSPEN(y0/(y0-y1)) - CSPEN((y0-DCMPLX(1D0))/(y0-y1))

        reseta = (0d0,0d0)
        if (dreal(y0).ne.1d0) then
     
        reseta = ETAD(-y1,1d0/(y0-y1))*cdlog(y0/(y0-y1)) 
     &         - ETAD(DCMPLX(1D0)-y1,1d0/(y0-y1))*cdlog((y0-DCMPLX(1D0))/(y0-y1))

        endif

        S2int = res + reseta

      endif

c if a,b are complex, one might need some additional eta functions (not correct below!)
c       zz = a*y0 + b
c
c       if (dreal(a).le.0d0 .and. aimag(b).ge.0d0 .and. aimag(zz).le.0d0) then
c             reseta2 = - 2*pi*cI
c       endif
c        
c       reseta = reseta + reseta2*cdlog((y0-DCMPLX(1D0))/y0)

      END

************************************************************************
        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
************************************************************************
        IMPLICIT REAL*8 (A-Z)
	REAL*8 M(4),P(4,4),K(4,4)
        real*8 pi,eps
        real*8 im1,im2
        real*8 m1,m2,m3,m4
        real*8 m02,m12,m22,m32,m42
        real*8 mm0,mm1,mm2,mm3,mm4 
        real*8 p1,p2,p3,p4,p12,p23
        real*8 q0,q1,q2,q3,q4,q00,q12,q23
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D04,CSPEN,ETAD,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	INTEGER I,J

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)
		
	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D04 = D04 + (-1D0)**(I+J)*(
     *		CSPEN(A1)+ETAD(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPEN(A2)+ETAD(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )

	   D04 = D04 
     &         + (3D0-2D0*I)*(
     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
     *		       + ETAD( RS(1,3),1D0/RS(2,4) )                  )
     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	     + ETAD( RS(1,3),1D0/RS(2,4) )
     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 
     &          + (3D0-2D0*I)*(
     *		     ETAD(-XX(I),1D0/R(2,4))
     *		      *( LOG(XX(I)/R(2,4)) + L1 )
     *		    +ETAD(-XX(I),R(1,3))
     *		      *( LOG(R(1,3)*XX(I)) + L2 )
     *		    -( ETAD(-XX(I),R(1,3)/R(2,4))
     *		      +ETAD(R(1,3),1D0/R(2,4)) )
     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	    +ETAD(R(1,3),1D0/R(2,4))
     *		      *ETAD(-XX(I),-R(1,3)/R(2,4))
     *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF
	
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       -CSPEN(1D0+SS(1)*X(I,1))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+X(I,1)/SS(1))
     *	       +ETAD(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETAD(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETAD(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETAD(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPEN(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPEN(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETAD(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETAD(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETAD(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
              pi = 4.D0*datan(1.D0)
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       +ETAD(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETAD(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPEN(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETAD(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETAD(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        integer k
        COMPLEX*16 CSPEN,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)
        data B/
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=54513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
c      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPEN=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPEN=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/dcmplx(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPEN=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/dcmplx(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/dcmplx(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/dcmplx(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPEN=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAD(C1,C2)                                            
***********************************************************************
*       COMPLEX ETAD-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETAD,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETAD = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETAD = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETAD = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETAD-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETAD,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETAD(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            

***********************************************************************
        FUNCTION KappaFun(A,B,C)                                            
***********************************************************************
*       for IR divergent box    				      *
*---------------------------------------------------------------------*
*       20.3.09  KK						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,KappaFun
	  
	if(dabs(dreal(A-(B-C)**2)).lt.10d-10) then
	  KappaFun = -1
	else
	  KappaFun = -(cdsqrt(dcmplx(1d0-4*B*C/(A-(B-C)**2)))-1d0)/(cdsqrt(dcmplx(1d0-4*B*C/(A-(B-C)**2)))+1d0)		
	endif

        END                                                            

!*******************************************************************
! Returns the solutions  x1,x2  to the equation  aa*x^2+bb*x+cc=0
! Also returns  dd = aa*(x1-x2)
! If  imode=/=0  it uses  dd  as input as value of  sqrt(b^2-4*a*c)
!*******************************************************************
      subroutine Qsol( x1,x2 ,dd ,aa,bb,cc,imode)
      
      implicit none
      complex*16 x1,x2 ,dd ,aa,bb,cc,czero
      complex*16 qq,hh
      double precision r1,r2
      integer imode
      
      czero = (0d0,0d0)
      
      if (aa.eq.czero) then
        if (bb.eq.czero) then
          x1= (0d0,0d0)
          x2= (0d0,0d0)
          dd= (0d0,0d0)
        else
          x1 = -cc/bb
          x2 = x1
          dd = bb                  
        endif
      else if (cc.eq.czero) then
        dd = -bb
        x1 = dd/aa
        x2 = (0d0,0d0)
      else
        if (imode.eq.0) dd = cdsqrt(bb*bb - 4*aa*cc)
        qq = -bb+dd
        hh = -bb-dd
        r1 = cdabs(qq)
        r2 = cdabs(hh)
        if (r1.ge.r2) then
          x1 = qq/(2*aa)
          x2 = (2*cc)/qq
        else
          qq = hh
          x2 = qq/(2*aa)
          x1 = (2*cc)/qq
        endif            
      end if
      end
      
***********************************************************************
        FUNCTION ETAD3(C1,C2,C3)                                            
***********************************************************************
*       COMPLEX ETAD-FUNKTION                                           
*---------------------------------------------------------------------*
*       23.05.13    K.Kovarik                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETAD3,C1,C2,C3                                           
        REAL*8     PI,IM1,IM2,IM12                                     

        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C3)                                          

        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETAD3 = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETAD3 = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETAD3 = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            
      
***********************************************************************
      subroutine setepsilon()
***********************************************************************
*       subroutine to set machine epsilon
*---------------------------------------------------------------------*
*       29.05.13    K.Kovarik                                       
***********************************************************************
      implicit none
      double precision MACHEPS

      MACHEPS = 1.D0
100   CONTINUE
      MACHEPS = MACHEPS / 2.D0
      IF ( 1.D0 + MACHEPS / 2.D0 .EQ. 1.D0 ) GOTO 110
      GO TO 100
110   CONTINUE

      END      

***********************************************************************
      function MACHEPS()
***********************************************************************
*       subroutine to set machine epsilon
*---------------------------------------------------------------------*
*       29.05.13    K.Kovarik                                       
***********************************************************************
      implicit none
      double precision MACHEPS

      MACHEPS = 1.D0
100   CONTINUE
      MACHEPS = MACHEPS / 2.D0
      IF ( 1.D0 + MACHEPS / 2.D0 .EQ. 1.D0 ) GOTO 110
      GO TO 100
110   CONTINUE

      END     
   
      
c***********************************************************************
c	Function LLi2(a,b)
c***********************************************************************
c*     analytical continued di-logarithm with 2 arguments Li_2 (1-ab)
c*     Eq. (2.14) in 1005.2076
c*---------------------------------------------------------------------*
c*       09.07.13    PS
c***********************************************************************
c	
c	implicit none
c	complex*16 a,b,LLi2, ETAD, cspen
c	
c	LLi2 = cspen(dcmplx(1d0-a*b)) + ETAD(a,b)*cdlog(dcmplx(1d0-a*b))
c	
c	end

***********************************************************************
	Function LLi2(z1,z2)
***********************************************************************
*     analytical continued di-logarithm with 2 arguments Li_2 (1-z1*z2)
*     Eq. (2.12) in 1005.2076
*---------------------------------------------------------------------*
*       09.07.13    KK                                       
***********************************************************************
	
	implicit none
	complex*16 z1,z2,z12,LLi2, ETAD, cspen
      double precision acc,eps

#include "../DMNLO_LoopIntegrals.h"

      acc = 1d-12
      eps = 1d-25

	z12 = z1*z2
	
      if( abs(z12) .le. 1d0 ) then
c Li2(1)            
	      if(abs(z12).lt.eps) then
	         LLi2 = Lpi**2
c Li2(0)            
	      else if (abs(z12 - 1d0).lt.acc) then
	         LLi2 = 0d0
c Li2(1-z1*z2) for |z12| < 1   
	      else
	         LLi2 = cspen(1d0-z12) 
     &         + (cdlog(z12) - cdlog(z1)  - cdlog(z2))*cdlog(1d0-z12)
	      endif
c Li2(1-z1*z2) for |z123| > 1   
	else 
	         LLi2 = cspen(1d0-z12) 
     &         + (cdlog(z12) - cdlog(z1) - cdlog(z2))*cdlog(1d0-z12)
     &         - (cdlog(-z12) - 0.5d0*(cdlog(z12)+cdlog(z1)+cdlog(z2)) )
	endif
		
	end


***********************************************************************
	Function LLi3(z1,z2,z3)
***********************************************************************
*     analytical continued di-logarithm with 3 arguments Li_2 (1-z1*z2*z3)
*     Eq. (2.12) in 1005.2076
*---------------------------------------------------------------------*
*       09.07.13    KK                                       
***********************************************************************
	
	implicit none
	complex*16 z1,z2,z3,z123,LLi3, ETAD, cspen
      double precision acc,eps

#include "../DMNLO_LoopIntegrals.h"

      acc = 1d-12
      eps = 1d-25

	z123 = z1*z2*z3
	
      if( abs(z123) .le. 1d0 ) then
c Li2(1)            
	      if(abs(z123).lt.eps) then
	         LLi3 = Lpi**2
c Li2(0)            
	      else if (abs(z123 - 1d0).lt.acc) then
	         LLi3 = 0d0
c Li2(1-z1*z2*z3) for |z123| < 1   
	      else
	         LLi3 = cspen(1d0-z123) 
     &         + (cdlog(z123) - cdlog(z1)  - cdlog(z2) - cdlog(z3))*cdlog(1d0-z123)
	      endif
c Li2(1-z1*z2*z3) for |z123| > 1   
	else 
	         LLi3 = cspen(1d0-z123) 
     &         + (cdlog(z123) - cdlog(z1) - cdlog(z2) - cdlog(z3))*cdlog(1d0-z123)
     &         - (cdlog(-z123) - 0.5d0*(cdlog(z123)+cdlog(z1)+cdlog(z2)+cdlog(z3)) )
	endif
		
	end

	
!==================================================================================================
!=========== Calculates C0, C1, C2 and C00 in the case of a vanishing Gram determinant ============
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine alternativeCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2)
	implicit none
	
!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "../DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	integer i, k
	double precision p1p1pat, p2p2pat, p1p2pat, m1m1, m2m2, m0m0, f1, f2
	double precision Mandelvar
	double precision x1(2), x2(2), Z11(2), Z12(2), Z21(2), Z22(2), DetZ(2), Y(2)
	complex *16 R7(2), R8(2)
	complex *16 B012, B112, C1_alt(2), C2_alt(2), C00_alt(2), Caa_alt(2)
	complex *16 C11_alt(2), C12_alt(2), C22_alt(2)
	complex *16 Ca_alt(2), B001, B002, B101, B102
	complex *16 B0001, B0002, B0012, B1101, B1102, B1112
				
!----------------------- Preparations -------------------------------------------------------------
	do i = 1,2
		C1_alt(i) = 0d0 
		C2_alt(i) = 0d0
		Ca_alt(i) = 0d0
		C00_alt(i) = 0d0
		Caa_alt(i) = 0d0
		C11_alt(i) = 0d0
		C12_alt(i) = 0d0
		C22_alt(i) = 0d0
	enddo
	
	call myBget(p1p1pat,m0m0,m1m1)
	B001 = xb0
	B101 = xb1
	B0001 = xb00
	B1101 = xb11
	
	call myBget(p2p2pat,m0m0,m2m2)
	B002 = xb0
	B102 = xb1
	B0002 = xb00
	B1102 = xb11
	
	call myBget(Mandelvar,m1m1,m2m2)
	B012 = xb0
	B112 = xb1
	B0012 = xb00
	B1112 = xb11
	
	f1 = p1p1pat - m1m1 + m0m0
	f2 = p2p2pat - m2m2 + m0m0
	p1p2pat = -0.5d0*(Mandelvar - p1p1pat - p2p2pat)
		
!----------------------- Select special case ------------------------------------------------------

! t-channel with both momenta equal and zero momentum transfer C(p**2, 0, p**2, m0**2, m1**2, m2**2)	
	if ((Mandelvar.eq.0d0).and.(p1p1pat.eq.p2p2pat)) then
				
		if(p1p1pat.eq.0d0) then
			Ca_alt(1) = 0d0
			Ca_alt(2) = 0d0
			Caa_alt(1) = 0d0
			Caa_alt(2) = 0d0
			if(EpsPole.ge.0d0) then
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0)
			endif
		else
			Ca_alt(1) = 1/(2d0*p1p1pat)*(B001-B012-f2*xc0)
			Ca_alt(2) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
			if(EpsPole.ge.0d0) then
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B102 - 2d0*f1*Ca_alt(1) - 0.5d0) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B101 - 2d0*f2*Ca_alt(2) - 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B102 + f1*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B101 + f2*Ca_alt(2) + 1d0)
			else
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B102 - 2d0*f1*Ca_alt(1)) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B101 - 2d0*f2*Ca_alt(2))
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B102 + f1*Ca_alt(1))
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B101 + f2*Ca_alt(2))
			endif
		endif
		
		xc1 = Ca_alt(1)/2d0
		xc2 = Ca_alt(1)/2d0
		xc11 = Caa_alt(1)/4d0
		xc12 = Caa_alt(1)/4d0
		xc22 = Caa_alt(1)/4d0
		xc00 = C00_alt(1)
		
! s-channel or u-channel with one momentum zero and zero momentum transfer C(p**2, p**2, 0, m0**2, m1**2, m2**2)
	elseif((Mandelvar.eq.p1p1pat).and.(p2p2pat.eq.0d0)) then
		
		if(p1p1pat.eq.0d0) then
			Ca_alt(1) = 0d0
			Ca_alt(2) = 0d0
			Caa_alt(1) = 0d0
			Caa_alt(2) = 0d0
			if(EpsPole.ge.0d0) then
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0)
			endif	
		else
			Ca_alt(1) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
			Ca_alt(2) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
			if(EpsPole.ge.0d0) then
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(1) - 0.5d0) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(2) - 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(2) + 1d0)
			else
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(1)) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(2))
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(1))
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(2))
			endif
		endif
			
		xc1 = Ca_alt(1)
		xc2 = 0d0
		xc11 = Caa_alt(1)
		xc12 = 0d0
		xc22 = 0d0
		xc00 = C00_alt(1)
			
! s-channel or u-channel with one momentum zero and zero momentum transfer C(0, p**2, p**2, m0**2, m1**2, m2**2)
	elseif((Mandelvar.eq.p2p2pat).and.(p1p1pat.eq.0d0)) then
		
		if(p2p2pat.eq.0d0) then
			Ca_alt(1) = 0d0
			Ca_alt(2) = 0d0
			Caa_alt(1) = 0d0
			Caa_alt(2) = 0d0
			if(EpsPole.ge.0d0) then
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0)
			endif
		else
			Ca_alt(1) = 1/(2d0*p2p2pat)*(B001-B012-f2*xc0)
			Ca_alt(2) = 1/(2d0*p2p2pat)*(B001-B012-f2*xc0)
			if(EpsPole.ge.0d0) then
				Caa_alt(1) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(1) + 0.5d0) 
				Caa_alt(2) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(2) + 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(2) + 1d0)
			else
				Caa_alt(1) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(1)) 
				Caa_alt(2) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(2))
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(1))
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(2))
			endif
		endif
			
		xc1 = 0d0
		xc2 = Ca_alt(1)
		xc11 = 0d0
		xc12 = 0d0
		xc22 = Caa_alt(1)
		xc00 = C00_alt(1)
				
! s-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)+m(p2))**2, p2**2, m0**2, m1**2, m2**2)
	elseif((Mandelvar.eq.((dsqrt(p1p1pat) + dsqrt(p2p2pat))**2).or.(abs(Mandelvar - ((dsqrt(p1p1pat) + dsqrt(p2p2pat))**2)).lt.0.01d0)).and.(p1p1pat.ne.0d0).and.(p2p2pat.ne.0d0)) then

		x1(1) = p2p2pat
		x1(2) = -p1p2pat
		x2(1) = -p1p2pat
		x2(2) = p1p1pat
				
		do k=1,2
		
			Y(k) = x1(k)*f1 + x2(k)*f2
			
			Z11(k) = Y(k) + x1(k)/2d0*f1
			Z12(k) = x1(k)/2d0*f2
			Z21(k) = x2(k)/2d0*f1
			Z22(k) = Y(k) + x2(k)/2d0*f2
			
			DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
			
			if(EpsPole.ge.0d0) then
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012 + 0.5d0)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012 + 0.5d0))
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k) + 1d0)
			else
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012))
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k))
			endif
			
		enddo
		
! Mathematica results
		if(EpsPole.ge.0d0) then
	
	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + p2p2pat/9. - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(p1p1pat/9. - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(p1p2pat/18. + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          p2p2pat/18. - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))
 
	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	5*f2**2*p1p1pat*p1p2pat**2 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	3*f1*f2*p1p2pat**3 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	5*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 3*f2**2*p1p1pat*p1p2pat*p2p2pat + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 3*f1**2*p1p2pat**2*p2p2pat - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 3*f1*f2*p1p2pat**2*p2p2pat - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 5*f1*f2*p1p1pat*p2p2pat**2 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 5*f1**2*p1p2pat*p2p2pat**2 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))
 
        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 + 
     -        	f1*f2*p1p1pat*p1p2pat**2 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -       	4*f2**2*p1p1pat*p1p2pat**2 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 + 
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 9*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 5*f1**2*p1p1pat*p2p2pat**2 + 
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 - 
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))
	

	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 2*p1p2pat - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat - 
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) - 
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(2*p1p2pat - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) - 
     -          (-p1p1pat - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + p1p2pat + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) - 
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))

        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	5*f2**2*p1p1pat**2*p1p2pat - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat + 
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 5*f1*f2*p1p1pat*p1p2pat**2 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 3*f2**2*p1p1pat*p1p2pat**2 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 - 
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 5*f1*f2*p1p2pat**3 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat + 
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 3*f1*f2*p1p1pat**2*p2p2pat + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat + 
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 3*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 3*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 5*f1**2*p1p2pat**2*p2p2pat + 
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat - 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat + 
     -        	f1*f2*p1p1pat**2*p1p2pat + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	4*f2**2*p1p1pat**2*p1p2pat + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat + 
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - f1**2*p1p1pat*p1p2pat**2 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 - 
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 9*f1*f2*p1p1pat*p1p2pat**2 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 5*f1**2*p1p2pat**3 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))
		
 		else
			
	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + 0d0 - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(0d0 - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(0d0 + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          0d0 - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))
 
	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	0d0 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	0d0 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	0d0 + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 0d0 + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 0d0 - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 0d0 - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 0d0 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 0d0 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))
	 	
        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 +
     -        	0d0 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 -
     -       	0d0 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 +
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - 0d0 -
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 0d0 -
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat -
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 0d0 +
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 -
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat +
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))

	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 0d0 - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat -
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) -
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(0d0 - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) -
     -          (0d0 - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + 0d0 + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) -
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))

        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat -
     -        	0d0 - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat +
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 0d0 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 +
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 0d0 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 -
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 0d0 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 +
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat +
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 0d0 + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat +
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 0d0 - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat -
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 0d0 -
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat -
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 0d0 +
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat -
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat +
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat +
     -        	0d0 + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat -
     -        	0d0 + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat +
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - 0d0 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 -
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 0d0 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 +
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 0d0 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 +
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat +
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat +
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))
			
		endif
		
		xc1 = C1_alt(1)
		xc2 = C2_alt(1)
		xc00 = C00_alt(1)
		xc11 = C11_alt(1)
		xc12 = C12_alt(1)
		xc22 = C22_alt(1)
		
! u-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)-m(p2))**2, p2**2, m0**2, m1**2, m2**2)
	elseif((Mandelvar.eq.((dsqrt(p1p1pat) - dsqrt(p2p2pat))**2).or.(abs(Mandelvar - ((dsqrt(p1p1pat) - dsqrt(p2p2pat))**2)).lt.0.01d0)).and.(p1p1pat.ne.0d0).and.(p2p2pat.ne.0d0)) then
		
		x1(1) = p2p2pat
		x1(2) = -p1p2pat
		x2(1) = -p1p2pat
		x2(2) = p1p1pat
				
		do k=1,2
		
			Y(k) = x1(k)*f1 + x2(k)*f2
			
			Z11(k) = Y(k) + x1(k)/2d0*f1
			Z12(k) = x1(k)/2d0*f2
			Z21(k) = x2(k)/2d0*f1
			Z22(k) = Y(k) + x2(k)/2d0*f2
			
			DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
			
			if(EpsPole.ge.0d0) then
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012 + 0.5d0)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012 + 0.5d0))
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k) + 1d0) 
			else
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012))
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k)) 
			endif
		enddo
		
! Mathematica results
		if(EpsPole.ge.0d0) then

	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + p2p2pat/9. - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(p1p1pat/9. - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(p1p2pat/18. + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          p2p2pat/18. - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))

	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	5*f2**2*p1p1pat*p1p2pat**2 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	3*f1*f2*p1p2pat**3 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	5*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 3*f2**2*p1p1pat*p1p2pat*p2p2pat + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 3*f1**2*p1p2pat**2*p2p2pat - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 3*f1*f2*p1p2pat**2*p2p2pat - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 5*f1*f2*p1p1pat*p2p2pat**2 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 5*f1**2*p1p2pat*p2p2pat**2 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))

        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 + 
     -        	f1*f2*p1p1pat*p1p2pat**2 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -       	4*f2**2*p1p1pat*p1p2pat**2 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 + 
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 9*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 5*f1**2*p1p1pat*p2p2pat**2 + 
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 - 
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))


	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 2*p1p2pat - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat - 
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) - 
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(2*p1p2pat - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) - 
     -          (-p1p1pat - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + p1p2pat + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) - 
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))

        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	5*f2**2*p1p1pat**2*p1p2pat - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat + 
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 5*f1*f2*p1p1pat*p1p2pat**2 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 3*f2**2*p1p1pat*p1p2pat**2 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 - 
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 5*f1*f2*p1p2pat**3 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat + 
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 3*f1*f2*p1p1pat**2*p2p2pat + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat + 
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 3*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 3*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 5*f1**2*p1p2pat**2*p2p2pat + 
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat - 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat + 
     -        	f1*f2*p1p1pat**2*p1p2pat + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	4*f2**2*p1p1pat**2*p1p2pat + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat + 
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - f1**2*p1p1pat*p1p2pat**2 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 - 
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 9*f1*f2*p1p1pat*p1p2pat**2 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 5*f1**2*p1p2pat**3 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

 		else
	
	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + 0d0 - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(0d0 - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(0d0 + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          0d0 - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))

	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	0d0 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	0d0 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	0d0 + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 0d0 + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 0d0 - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 0d0 - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 0d0 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 0d0 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))

        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 +
     -        	0d0 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 -
     -       	0d0 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 +
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - 0d0 -
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 0d0 -
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat -
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 0d0 +
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 -
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat +
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))

	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 0d0 - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat -
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) -
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(0d0 - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) -
     -          (0d0 - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + 0d0 + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) -
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))

        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat -
     -        	0d0 - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat +
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 0d0 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 +
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 0d0 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 -
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 0d0 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 +
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat +
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 0d0 + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat +
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 0d0 - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat -
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 0d0 -
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat -
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 0d0 +
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat -
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat +
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat +
     -        	0d0 + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat -
     -        	0d0 + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat +
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - 0d0 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 -
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 0d0 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 +
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 0d0 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 +
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat +
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat +
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))
	
		endif

		xc1 = C1_alt(1)
		xc2 = C2_alt(1)
		xc00 = C00_alt(1)
		xc11 = C11_alt(1)
		xc12 = C12_alt(1)
		xc22 = C22_alt(1)
		
	else
	
		write(*,*)'Error in alternative C tensor determination - Case not found'

	endif
				
!======================= END OF PROGRAM ===========================================================

	end
	

!==================================================================================================
!============= Calculates D functions in the case of a vanishing Gram determinant =================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine alternativeDget(p10p10pat, p21p21pat, p32p32pat, p30p30pat, p20p20pat, p31p31pat, m0m0, m1m1, m2m2, m3m3)
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "../DMNLO_LoopIntegrals.h"
#include "../DMNLO_Model.h"

!----------------------- Declare variables --------------------------------------------------------
	double precision p10p10pat, p20p20pat, p30p30pat, p21p21pat, p32p32pat, p31p31pat
	double precision p1p2pat, p1p3pat, p2p3pat
	double precision m0m0, m1m1, m2m2, m3m3
	double precision f1, f2, f3
	double precision x1(2), x2(2), Z11(2), Z12(2), Z21(2), Z22(2), DetZ(2), Y(2)
	complex*16 W1(2), W2(2)
	complex*16 D1_alt(2), D2_alt(2), D3_alt(2), D00_alt(2), Da_alt(2)
	complex*16 D11_alt(2), D12_alt(2), D13_alt(2), D22_alt(2), D23_alt(2), D33_alt(2)
	complex*16 D3a_alt(2), Daa_alt(2)
	complex*16 Da_test(2), D3_test(2)
	complex*16 C0012, C0013, C0023, C0123
	complex*16 C1012, C1013, C1023, C1123
	complex*16 C2012, C2013, C2023, C2123
	complex*16 C00012, C00013, C00023, C00123
	complex*16 C11012, C11013, C11023, C11123
	complex*16 C12012, C12013, C12023, C12123
	complex*16 C22012, C22013, C22023, C22123
	integer i, k
	
c	open (unit=30,file="results_Noise.txt",action="write",status="old", position="append")
		
!----------------------- Preparations -------------------------------------------------------------
	do i = 1,2
		D1_alt(i) = 0d0
		D2_alt(i) = 0d0
		D3_alt(i) = 0d0
		Da_alt(i) = 0d0
		D00_alt(i) = 0d0
		D11_alt(i) = 0d0
		D12_alt(i) = 0d0
		D13_alt(i) = 0d0
		D22_alt(i) = 0d0
		D23_alt(i) = 0d0
		D33_alt(i) = 0d0
		Daa_alt(i) = 0d0
		D3a_alt(i) = 0d0
	enddo
			
	f1 = p10p10pat - m1m1 + m0m0
	f2 = p20p20pat - m2m2 + m0m0
	f3 = p30p30pat - m3m3 + m0m0

	p1p2pat = - 0.5d0*(p21p21pat - p10p10pat - p20p20pat)
	p1p3pat = - 0.5d0*(p31p31pat - p10p10pat - p30p30pat)
	p2p3pat = - 0.5d0*(p32p32pat - p20p20pat - p30p30pat)	

	call myCget(p10p10pat,p21p21pat, p20p20pat, m0m0, m1m1, m2m2)
	C0012 = xc0
	C1012 = xc1
	C2012 = xc2
	C00012 = xc00
	C11012 = xc11
	C12012 = xc12
	C22012 = xc22
	
	call myCget(p10p10pat,p31p31pat, p30p30pat, m0m0, m1m1, m3m3)
	C0013 = xc0
	C1013 = xc1
	C2013 = xc2
	C00013 = xc00
	C11013 = xc11
	C12013 = xc12
	C22013 = xc22
	
	call myCget(p20p20pat,p32p32pat, p30p30pat, m0m0, m2m2, m3m3)
	C0023 = xc0
	C1023 = xc1
	C2023 = xc2
	C00023 = xc00
	C11023 = xc11
	C12023 = xc12
	C22023 = xc22

	call myCget(p21p21pat,p32p32pat, p31p31pat, m1m1, m2m2, m3m3)
	C0123 = xc0
	C1123 = xc1
	C2123 = xc2
	C00123 = xc00
	C11123 = xc11
	C12123 = xc12
	C22123 = xc22

! Box B2 or B4 with zero momentum transfer D(p, p, p3, m0**2, m1**2, m2**2, m3**2) = (p**2, 0, (p3-p)**2, p3**2, p**2, (p3-p)**2, m0**2, m1**2, m2**2, m3**2)

	if ((p10p10pat.eq.p20p20pat).and.(p21p21pat.eq.0d0).and.(p31p31pat.eq.p32p32pat).and.(p10p10pat.ne.p30p30pat)) then
					
		if(p10p10pat.eq.0d0) then
				
			do i = 1,2
				D1_alt(i) = 0d0
				D2_alt(i) = 0d0
				D3_alt(i) = 1d0/(2d0*p30p30pat)*(C0012 - C0123 - f3*xd0)
				D11_alt(i) = 0d0
				D12_alt(i) = 0d0
				D13_alt(i) = 0d0
				D22_alt(i) = 0d0
				D23_alt(i) = 0d0
			enddo
		
			D00_alt(1) = 1d0/3d0*(C0123 + m0m0*xd0 - 0.5d0*(C2012 - C2123 - f3*D3_alt(1)))
			D00_alt(2) = 1d0/3d0*(C0023 + m1m1*xd0 - 0.5d0*(C2012 - C2123 - f3*D3_alt(2)))
			D33_alt(1) = 1d0/p30p30pat*(2d0/3d0*(C2012 -C2123 - f3*D3_alt(1)) - 1d0/3d0*(C0123 + m0m0*xd0))
			D33_alt(2) = 1d0/p30p30pat*(2d0/3d0*(C2012 -C2123 - f3*D3_alt(2)) - 1d0/3d0*(C0023 + m1m1*xd0))
			
		else
		
			x1(1) = p10p10pat
			x1(2) = - p1p3pat
			x2(1) = - p1p3pat
			x2(2) = p30p30pat
						
			do k = 1,2
		
				Y(k) = x1(k)*f3 + x2(k)*f1
	
				Z11(k) = Y(k) + x1(k)/2d0*f3
				Z12(k) = x1(k)/2d0*f1
				Z21(k) = x2(k)/2d0*f3
				Z22(k) = Y(k) + x2(k)/2d0*f1
								
				DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
							
				W1(k) = x1(k)*(-C2123 - m0m0*xd0 + 0.5*(C1023 - C0123)) + x2(k)*(C2023 - C2123)
				W2(k) = x1(k)*(C1012 + C2012 + C0123 + C2123) + x2(k)*(3d0/2d0*C1023 + C2123 - m0m0*xd0 + 0.5d0*C0123)
				
				D3_alt(k) = 1d0/DetZ(k)*(Z22(k)*W1(k) - Z12(k)*W2(k))
				Da_alt(k) = 1d0/DetZ(k)*(-Z21(k)*W1(k) + Z11(k)*W2(k))
				
c				write(*,*)'DetZ',DetZ(k)
c				write(*,*)'-Z21(k)*W1(k)',-Z21(k)*W1(k)
c				write(*,*)'Z11(k)*W2(k)',Z11(k)*W2(k)
c				write(*,*)'-Z21(k)*W1(k) + Z11(k)*W2(k)',-Z21(k)*W1(k) + Z11(k)*W2(k)				
c				write(*,*)'Da_alt(k)', Da_alt(k)
			
				D1_alt(k) = Da_alt(k)/2d0
				D2_alt(k) = Da_alt(k)/2d0
				D00_alt(k) = 0.5d0*(m0m0*xd0 + 0.5d0*(f3*D3_alt(k) + f1*Da_alt(k) - C1023 + C0123))
			
			enddo
			
c	        write (30,*),REAL(MUE),REAL(DetZ(1)),REAL(-Z21(1)*W1(1)),REAL(Z11(1)*W2(1)),REAL(-Z21(1)*W1(1) + Z11(1)*W2(1)),REAL(Da_alt(1))		
		
! Mathematica results for D3 and Da				
C 			D3_test(1) = -((-((-2d0*f3*p1p3pat + 3d0*f1*p30p30pat)*(-(C0123*p10p10pat)/4d0 + (C1023*p10p10pat)/4d0 - (C2123*p10p10pat)/2d0 - (C2023*p1p3pat)/2d0 +
C      -     	(C2123*p1p3pat)/2d0 - (m0m0*p10p10pat*xd0)/2d0)) -
C      -        (f1*p10p10pat*(2d0*C0123*p1p3pat + 2d0*C1012*p1p3pat + 2d0*C2012*p1p3pat + 2d0*C2123*p1p3pat - C0123*p30p30pat - 3d0*C1023*p30p30pat - 2d0*C2123*p30p30pat +
C      -        2d0*m0m0*p30p30pat*xd0))/4d0)/((-3d0*f3**2*p10p10pat*p1p3pat)/2d0 + f1*f3*p1p3pat**2 + 2d0*f1*f3*p10p10pat*p30p30pat - (3d0*f1**2*p1p3pat*p30p30pat)/2d0))
C
C 			Da_test(1) = -((-3d0*C0123*f3*p10p10pat*p1p3pat - 3d0*C1012*f3*p10p10pat*p1p3pat - 3d0*C2012*f3*p10p10pat*p1p3pat - 3d0*C2123*f3*p10p10pat*p1p3pat + 2d0*C0123*f1*p1p3pat**2 +
C      -        2d0*C1012*f1*p1p3pat**2 + 2d0*C2012*f1*p1p3pat**2 + 2d0*C2123*f1*p1p3pat**2 + 2d0*C0123*f3*p10p10pat*p30p30pat + 4d0*C1023*f3*p10p10pat*p30p30pat +
C      -        4d0*C2123*f3*p10p10pat*p30p30pat - C0123*f1*p1p3pat*p30p30pat - 3d0*C1023*f1*p1p3pat*p30p30pat - 2d0*C2123*f1*p1p3pat*p30p30pat + C2023*f3*p1p3pat*p30p30pat -
C      -        C2123*f3*p1p3pat*p30p30pat - 2d0*f3*m0m0*p10p10pat*p30p30pat*xd0 + 2d0*f1*m0m0*p1p3pat*p30p30pat*xd0)/
C      -		(3d0*f3**2*p10p10pat*p1p3pat - 2d0*f1*f3*p1p3pat**2 - 4d0*f1*f3*p10p10pat*p30p30pat + 3d0*f1**2*p1p3pat*p30p30pat))
C
C 			D3_test(2) = -((-(f1*p1p3pat*(-2d0*C0123*p10p10pat - 2d0*C1012*p10p10pat - 2d0*C2012*p10p10pat - 2d0*C2123*p10p10pat + C0123*p1p3pat + 3d0*C1023*p1p3pat +
C      -        2d0*C2123*p1p3pat - 2d0*m0m0*p1p3pat*xd0)) - (2d0*f3*p10p10pat - 3d0*f1*p1p3pat)*
C      -        (-(C0123*p1p3pat) + C1023*p1p3pat - 2d0*C2123*p1p3pat - 2d0*C2023*p30p30pat + 2d0*C2123*p30p30pat - 2d0*m0m0*p1p3pat*xd0))/
C      -		(f1*f3*p1p3pat**2 - (2d0*f3*p10p10pat - 3d0*f1*p1p3pat)*(-3d0*f3*p1p3pat + 2d0*f1*p30p30pat)))
C
C 			Da_test(2) = -((-3d0*C0123*f3*p10p10pat*p1p3pat - 3d0*C1012*f3*p10p10pat*p1p3pat - 3d0*C2012*f3*p10p10pat*p1p3pat - 3d0*C2123*f3*p10p10pat*p1p3pat + 2d0*C0123*f3*p1p3pat**2 +
C      -        4d0*C1023*f3*p1p3pat**2 + 4d0*C2123*f3*p1p3pat**2 + 2d0*C0123*f1*p10p10pat*p30p30pat + 2d0*C1012*f1*p10p10pat*p30p30pat + 2d0*C2012*f1*p10p10pat*p30p30pat +
C      -        2d0*C2123*f1*p10p10pat*p30p30pat - C0123*f1*p1p3pat*p30p30pat - 3d0*C1023*f1*p1p3pat*p30p30pat - 2d0*C2123*f1*p1p3pat*p30p30pat + C2023*f3*p1p3pat*p30p30pat -
C      -        C2123*f3*p1p3pat*p30p30pat - 2d0*f3*m0m0*p1p3pat**2*xd0 + 2d0*f1*m0m0*p1p3pat*p30p30pat*xd0)/
C      -		(3d0*f3**2*p10p10pat*p1p3pat - 4d0*f1*f3*p1p3pat**2 - 2d0*f1*f3*p10p10pat*p30p30pat + 3d0*f1**2*p1p3pat*p30p30pat))
				
! Mathematica results for D33, D3a and Daa						
			D33_alt(1) = -((-(((C12023*p10p10pat)/3. - (C2123*p10p10pat)/3. - (C22123*p10p10pat)/2. - (2*D3_alt(1)*m0m0*p10p10pat)/3. - (C22023*p1p3pat)/2. + 
     -          (C22123*p1p3pat)/2.)*((f3**2*p1p3pat**2)/3. - (5*f1*f3*p1p3pat*p30p30pat)/6. + (5*f1**2*p30p30pat**2)/9.)) - 
     -         	(f1*p10p10pat*(-(((f3*p1p3pat)/2. - (5*f1*p30p30pat)/6.)*
     -          ((C0123*p1p3pat)/6. + (C11023*p1p3pat)/6. + (2*C2123*p1p3pat)/3. + (C22123*p1p3pat)/2. - (Da_alt(1)*m0m0*p1p3pat)/3. - (2*C12023*p30p30pat)/3. - 
     -          (C2123*p30p30pat)/3. - (C22123*p30p30pat)/2. + (D3_alt(1)*m0m0*p30p30pat)/3.)) - 
     -          (f1*p1p3pat*((C0123*p1p3pat)/2. - (C11012*p1p3pat)/2. - C12012*p1p3pat + C2123*p1p3pat - (C22012*p1p3pat)/2. + (C22123*p1p3pat)/2. - 
     -          (C0123*p30p30pat)/6. + (5*C11023*p30p30pat)/6. - (2*C2123*p30p30pat)/3. - (C22123*p30p30pat)/2. - (2*Da_alt(1)*m0m0*p30p30pat)/3.))/6.))/3.)/
     -       	(-(((-5*f3*p10p10pat)/6. + (f1*p1p3pat)/2.)*((f3**2*p1p3pat**2)/3. - (5*f1*f3*p1p3pat*p30p30pat)/6. + (5*f1**2*p30p30pat**2)/9.)) - 
     -         	(f1*p10p10pat*(-(f3**2*p1p3pat*p30p30pat)/12. + (5*f1*f3*p30p30pat**2)/36.))/3.))
 
     		D3a_alt(1) = -(5*C0123*f1*f3*p10p10pat*p1p3pat**2 - 5*C11012*f1*f3*p10p10pat*p1p3pat**2 - 10*C12012*f1*f3*p10p10pat*p1p3pat**2 + 
     -        	10*C2123*f1*f3*p10p10pat*p1p3pat**2 - 5*C22012*f1*f3*p10p10pat*p1p3pat**2 + 5*C22123*f1*f3*p10p10pat*p1p3pat**2 + 5*C0123*f3**2*p10p10pat*p1p3pat**2 + 
     -        	5*C11023*f3**2*p10p10pat*p1p3pat**2 + 20*C2123*f3**2*p10p10pat*p1p3pat**2 + 15*C22123*f3**2*p10p10pat*p1p3pat**2 - 10*Da_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat**2 - 
     -        	3*C0123*f1**2*p1p3pat**3 + 3*C11012*f1**2*p1p3pat**3 + 6*C12012*f1**2*p1p3pat**3 - 6*C2123*f1**2*p1p3pat**3 + 3*C22012*f1**2*p1p3pat**3 - 
     -        	3*C22123*f1**2*p1p3pat**3 - 3*C0123*f1*f3*p1p3pat**3 - 3*C11023*f1*f3*p1p3pat**3 - 12*C2123*f1*f3*p1p3pat**3 - 9*C22123*f1*f3*p1p3pat**3 + 
     -        	6*Da_alt(1)*f1*f3*m0m0*p1p3pat**3 - 10*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 40*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 
     -        	30*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 18*C12023*f3**2*p10p10pat*p1p3pat*p30p30pat - 12*C2123*f3**2*p10p10pat*p1p3pat*p30p30pat - 
     -        	18*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat + 10*Da_alt(1)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat + 6*D3_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat*p30p30pat + 
     -        	6*C0123*f1**2*p1p3pat**2*p30p30pat + 24*C2123*f1**2*p1p3pat**2*p30p30pat + 18*C22123*f1**2*p1p3pat**2*p30p30pat + 12*C12023*f1*f3*p1p3pat**2*p30p30pat + 
     -        	6*C2123*f1*f3*p1p3pat**2*p30p30pat + 9*C22123*f1*f3*p1p3pat**2*p30p30pat - 3*C22023*f3**2*p1p3pat**2*p30p30pat + 3*C22123*f3**2*p1p3pat**2*p30p30pat - 
     -        	6*Da_alt(1)*f1**2*m0m0*p1p3pat**2*p30p30pat - 6*D3_alt(1)*f1*f3*m0m0*p1p3pat**2*p30p30pat + 30*C12023*f1*f3*p10p10pat*p30p30pat**2 + 
     -        	20*C2123*f1*f3*p10p10pat*p30p30pat**2 + 30*C22123*f1*f3*p10p10pat*p30p30pat**2 - 10*D3_alt(1)*f1*f3*m0m0*p10p10pat*p30p30pat**2 - 
     -        	20*C12023*f1**2*p1p3pat*p30p30pat**2 - 10*C2123*f1**2*p1p3pat*p30p30pat**2 - 15*C22123*f1**2*p1p3pat*p30p30pat**2 + 5*C22023*f1*f3*p1p3pat*p30p30pat**2 - 
     -       	5*C22123*f1*f3*p1p3pat*p30p30pat**2 + 10*D3_alt(1)*f1**2*m0m0*p1p3pat*p30p30pat**2)/
     -     		(2.*(-10*f3**3*p10p10pat*p1p3pat**2 + 6*f1*f3**2*p1p3pat**3 + 24*f1*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*f1**2*f3*p1p3pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p30p30pat**2 + 10*f1**3*p1p3pat*p30p30pat**2))
 
 			Daa_alt(1) = -((-10*C0123*f3**2*p10p10pat*p1p3pat**2 + 10*C11012*f3**2*p10p10pat*p1p3pat**2 + 20*C12012*f3**2*p10p10pat*p1p3pat**2 - 
     -         	20*C2123*f3**2*p10p10pat*p1p3pat**2 + 10*C22012*f3**2*p10p10pat*p1p3pat**2 - 10*C22123*f3**2*p10p10pat*p1p3pat**2 + 6*C0123*f1*f3*p1p3pat**3 - 
     -         	6*C11012*f1*f3*p1p3pat**3 - 12*C12012*f1*f3*p1p3pat**3 + 12*C2123*f1*f3*p1p3pat**3 - 6*C22012*f1*f3*p1p3pat**3 + 6*C22123*f1*f3*p1p3pat**3 + 
     -         	9*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C11012*f1*f3*p10p10pat*p1p3pat*p30p30pat - 18*C12012*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -         	18*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C22012*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -         	5*C0123*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*C11023*f3**2*p10p10pat*p1p3pat*p30p30pat + 20*C2123*f3**2*p10p10pat*p1p3pat*p30p30pat + 
     -         	15*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat + 10*Da_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat*p30p30pat - 6*C0123*f1**2*p1p3pat**2*p30p30pat + 
     -         	6*C11012*f1**2*p1p3pat**2*p30p30pat + 12*C12012*f1**2*p1p3pat**2*p30p30pat - 12*C2123*f1**2*p1p3pat**2*p30p30pat + 6*C22012*f1**2*p1p3pat**2*p30p30pat - 
     -         	6*C22123*f1**2*p1p3pat**2*p30p30pat - 3*C0123*f1*f3*p1p3pat**2*p30p30pat + 9*C11023*f1*f3*p1p3pat**2*p30p30pat - 12*C2123*f1*f3*p1p3pat**2*p30p30pat - 
     -         	9*C22123*f1*f3*p1p3pat**2*p30p30pat - 6*Da_alt(1)*f1*f3*m0m0*p1p3pat**2*p30p30pat - 3*C0123*f1*f3*p10p10pat*p30p30pat**2 + 
     -         	15*C11023*f1*f3*p10p10pat*p30p30pat**2 - 12*C2123*f1*f3*p10p10pat*p30p30pat**2 - 9*C22123*f1*f3*p10p10pat*p30p30pat**2 - 
     -         	6*C12023*f3**2*p10p10pat*p30p30pat**2 - 4*C2123*f3**2*p10p10pat*p30p30pat**2 - 6*C22123*f3**2*p10p10pat*p30p30pat**2 - 
     -         	12*Da_alt(1)*f1*f3*m0m0*p10p10pat*p30p30pat**2 + 2*D3_alt(1)*f3**2*m0m0*p10p10pat*p30p30pat**2 + 2*C0123*f1**2*p1p3pat*p30p30pat**2 - 
     -         	10*C11023*f1**2*p1p3pat*p30p30pat**2 + 8*C2123*f1**2*p1p3pat*p30p30pat**2 + 6*C22123*f1**2*p1p3pat*p30p30pat**2 + 4*C12023*f1*f3*p1p3pat*p30p30pat**2 + 
     -         	2*C2123*f1*f3*p1p3pat*p30p30pat**2 + 3*C22123*f1*f3*p1p3pat*p30p30pat**2 - C22023*f3**2*p1p3pat*p30p30pat**2 + C22123*f3**2*p1p3pat*p30p30pat**2 + 
     -         	8*Da_alt(1)*f1**2*m0m0*p1p3pat*p30p30pat**2 - 2*D3_alt(1)*f1*f3*m0m0*p1p3pat*p30p30pat**2)/
     -       	(-10*f3**3*p10p10pat*p1p3pat**2 + 6*f1*f3**2*p1p3pat**3 + 24*f1*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*f1**2*f3*p1p3pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p30p30pat**2 + 10*f1**3*p1p3pat*p30p30pat**2))

	 		D33_alt(2) = -((-2*f1*p1p3pat*(-((3*f3*p10p10pat - 5*f1*p1p3pat)*
     -          (-(C0123*p10p10pat) - C11023*p10p10pat - 4*C2123*p10p10pat - 3*C22123*p10p10pat + 2*Da_alt(2)*m0m0*p10p10pat + 4*C12023*p1p3pat + 2*C2123*p1p3pat + 
     -          3*C22123*p1p3pat - 2*D3_alt(2)*m0m0*p1p3pat)) + f1*p10p10pat*
     -          (3*C0123*p10p10pat - 3*C11012*p10p10pat - 6*C12012*p10p10pat + 6*C2123*p10p10pat - 3*C22012*p10p10pat + 3*C22123*p10p10pat - C0123*p1p3pat + 
     -          5*C11023*p1p3pat - 4*C2123*p1p3pat - 3*C22123*p1p3pat - 4*Da_alt(2)*m0m0*p1p3pat)) - 
     -         	(-12*f3**2*p10p10pat**2 + 30*f1*f3*p10p10pat*p1p3pat - 20*f1**2*p1p3pat**2)*
     -          (2*C12023*p1p3pat - 2*C2123*p1p3pat - 3*C22123*p1p3pat - 4*D3_alt(2)*m0m0*p1p3pat - 3*C22023*p30p30pat + 3*C22123*p30p30pat))/
     -       	(-2*f1*f3*p1p3pat**2*(3*f3*p10p10pat - 5*f1*p1p3pat) - 
     -         	(-12*f3**2*p10p10pat**2 + 30*f1*f3*p10p10pat*p1p3pat - 20*f1**2*p1p3pat**2)*(-5*f3*p1p3pat + 3*f1*p30p30pat)))
 
     		D3a_alt(2) = -(5*C0123*f1*f3*p10p10pat**2*p1p3pat - 5*C11012*f1*f3*p10p10pat**2*p1p3pat - 10*C12012*f1*f3*p10p10pat**2*p1p3pat + 
     -        	10*C2123*f1*f3*p10p10pat**2*p1p3pat - 5*C22012*f1*f3*p10p10pat**2*p1p3pat + 5*C22123*f1*f3*p10p10pat**2*p1p3pat + 5*C0123*f3**2*p10p10pat**2*p1p3pat + 
     -        	5*C11023*f3**2*p10p10pat**2*p1p3pat + 20*C2123*f3**2*p10p10pat**2*p1p3pat + 15*C22123*f3**2*p10p10pat**2*p1p3pat - 10*Da_alt(2)*f3**2*m0m0*p10p10pat**2*p1p3pat - 
     -        	10*C0123*f1*f3*p10p10pat*p1p3pat**2 - 40*C2123*f1*f3*p10p10pat*p1p3pat**2 - 30*C22123*f1*f3*p10p10pat*p1p3pat**2 - 18*C12023*f3**2*p10p10pat*p1p3pat**2 - 
     -        	12*C2123*f3**2*p10p10pat*p1p3pat**2 - 18*C22123*f3**2*p10p10pat*p1p3pat**2 + 10*Da_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat**2 + 
     -        	6*D3_alt(2)*f3**2*m0m0*p10p10pat*p1p3pat**2 + 30*C12023*f1*f3*p1p3pat**3 + 20*C2123*f1*f3*p1p3pat**3 + 30*C22123*f1*f3*p1p3pat**3 - 10*D3_alt(2)*f1*f3*m0m0*p1p3pat**3 - 
     -        	3*C0123*f1**2*p10p10pat**2*p30p30pat + 3*C11012*f1**2*p10p10pat**2*p30p30pat + 6*C12012*f1**2*p10p10pat**2*p30p30pat - 
     -        	6*C2123*f1**2*p10p10pat**2*p30p30pat + 3*C22012*f1**2*p10p10pat**2*p30p30pat - 3*C22123*f1**2*p10p10pat**2*p30p30pat - 
     -        	3*C0123*f1*f3*p10p10pat**2*p30p30pat - 3*C11023*f1*f3*p10p10pat**2*p30p30pat - 12*C2123*f1*f3*p10p10pat**2*p30p30pat - 
     -        	9*C22123*f1*f3*p10p10pat**2*p30p30pat + 6*Da_alt(2)*f1*f3*m0m0*p10p10pat**2*p30p30pat + 6*C0123*f1**2*p10p10pat*p1p3pat*p30p30pat + 
     -        	24*C2123*f1**2*p10p10pat*p1p3pat*p30p30pat + 18*C22123*f1**2*p10p10pat*p1p3pat*p30p30pat + 12*C12023*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -        	6*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 3*C22023*f3**2*p10p10pat*p1p3pat*p30p30pat + 
     -        	3*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat - 6*Da_alt(2)*f1**2*m0m0*p10p10pat*p1p3pat*p30p30pat - 6*D3_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat - 
     -        	20*C12023*f1**2*p1p3pat**2*p30p30pat - 10*C2123*f1**2*p1p3pat**2*p30p30pat - 15*C22123*f1**2*p1p3pat**2*p30p30pat + 5*C22023*f1*f3*p1p3pat**2*p30p30pat - 
     -        	5*C22123*f1*f3*p1p3pat**2*p30p30pat + 10*D3_alt(2)*f1**2*m0m0*p1p3pat**2*p30p30pat)/
     -     		(2.*(-10*f3**3*p10p10pat**2*p1p3pat + 24*f1*f3**2*p10p10pat*p1p3pat**2 - 15*f1**2*f3*p1p3pat**3 + 6*f1*f3**2*p10p10pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p1p3pat*p30p30pat + 10*f1**3*p1p3pat**2*p30p30pat))
 
     		Daa_alt(2) = -((-10*C0123*f3**2*p10p10pat**2*p1p3pat + 10*C11012*f3**2*p10p10pat**2*p1p3pat + 20*C12012*f3**2*p10p10pat**2*p1p3pat - 
     -         	20*C2123*f3**2*p10p10pat**2*p1p3pat + 10*C22012*f3**2*p10p10pat**2*p1p3pat - 10*C22123*f3**2*p10p10pat**2*p1p3pat + 9*C0123*f1*f3*p10p10pat*p1p3pat**2 - 
     -         	9*C11012*f1*f3*p10p10pat*p1p3pat**2 - 18*C12012*f1*f3*p10p10pat*p1p3pat**2 + 18*C2123*f1*f3*p10p10pat*p1p3pat**2 - 9*C22012*f1*f3*p10p10pat*p1p3pat**2 + 
     -         	9*C22123*f1*f3*p10p10pat*p1p3pat**2 + 5*C0123*f3**2*p10p10pat*p1p3pat**2 - 15*C11023*f3**2*p10p10pat*p1p3pat**2 + 20*C2123*f3**2*p10p10pat*p1p3pat**2 + 
     -         	15*C22123*f3**2*p10p10pat*p1p3pat**2 + 10*Da_alt(2)*f3**2*m0m0*p10p10pat*p1p3pat**2 - 3*C0123*f1*f3*p1p3pat**3 + 15*C11023*f1*f3*p1p3pat**3 - 
     -         	12*C2123*f1*f3*p1p3pat**3 - 9*C22123*f1*f3*p1p3pat**3 - 6*C12023*f3**2*p1p3pat**3 - 4*C2123*f3**2*p1p3pat**3 - 6*C22123*f3**2*p1p3pat**3 - 
     -         	12*Da_alt(2)*f1*f3*m0m0*p1p3pat**3 + 2*D3_alt(2)*f3**2*m0m0*p1p3pat**3 + 6*C0123*f1*f3*p10p10pat**2*p30p30pat - 6*C11012*f1*f3*p10p10pat**2*p30p30pat - 
     -         	12*C12012*f1*f3*p10p10pat**2*p30p30pat + 12*C2123*f1*f3*p10p10pat**2*p30p30pat - 6*C22012*f1*f3*p10p10pat**2*p30p30pat + 
     -         	6*C22123*f1*f3*p10p10pat**2*p30p30pat - 6*C0123*f1**2*p10p10pat*p1p3pat*p30p30pat + 6*C11012*f1**2*p10p10pat*p1p3pat*p30p30pat + 
     -         	12*C12012*f1**2*p10p10pat*p1p3pat*p30p30pat - 12*C2123*f1**2*p10p10pat*p1p3pat*p30p30pat + 6*C22012*f1**2*p10p10pat*p1p3pat*p30p30pat - 
     -         	6*C22123*f1**2*p10p10pat*p1p3pat*p30p30pat - 3*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C11023*f1*f3*p10p10pat*p1p3pat*p30p30pat - 
     -         	12*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 6*Da_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat + 
     -         	2*C0123*f1**2*p1p3pat**2*p30p30pat - 10*C11023*f1**2*p1p3pat**2*p30p30pat + 8*C2123*f1**2*p1p3pat**2*p30p30pat + 6*C22123*f1**2*p1p3pat**2*p30p30pat + 
     -         	4*C12023*f1*f3*p1p3pat**2*p30p30pat + 2*C2123*f1*f3*p1p3pat**2*p30p30pat + 3*C22123*f1*f3*p1p3pat**2*p30p30pat - C22023*f3**2*p1p3pat**2*p30p30pat + 
     -         	C22123*f3**2*p1p3pat**2*p30p30pat + 8*Da_alt(2)*f1**2*m0m0*p1p3pat**2*p30p30pat - 2*D3_alt(2)*f1*f3*m0m0*p1p3pat**2*p30p30pat)/
     -       	(-10*f3**3*p10p10pat**2*p1p3pat + 24*f1*f3**2*p10p10pat*p1p3pat**2 - 15*f1**2*f3*p1p3pat**3 + 6*f1*f3**2*p10p10pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p1p3pat*p30p30pat + 10*f1**3*p1p3pat**2*p30p30pat))
		
			do k = 1,2
			
				D11_alt(k) = Daa_alt(k)/4d0
				D12_alt(k) = Daa_alt(k)/4d0
				D22_alt(k) = Daa_alt(k)/4d0
				D13_alt(k) = D3a_alt(k)/2d0
				D23_alt(k) = D3a_alt(k)/2d0
			
			enddo
			
		endif
		
		xd1 = D1_alt(1)
		xd2 = D2_alt(1)
		xd3 = D3_alt(1)
		xd00 = D00_alt(1)
		xd11 = D11_alt(1)
		xd12 = D12_alt(1)
		xd13 = D13_alt(1)
		xd22 = D22_alt(1)
		xd23 = D23_alt(1)
		xd33 = D33_alt(1)
		
	else

		write(*,*)'Error in alternative D tensor determination - Case not found'
		
	endif
	
c	close (unit=30)
	
!======================= END OF PROGRAM ===========================================================

	end

	
	 
