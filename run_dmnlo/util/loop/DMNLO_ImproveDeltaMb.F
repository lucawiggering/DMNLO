! ======================================================================================================
! This subroutine improves the SQCD resummation for the Higgs-b-bbar Coupling to NNLO. The most relevant 
! papers are arXiv:1001.1935, DOI: 10.1103/PhysRevD.68.115001 and DOI: 10.1016/S0550-3213(00)00146-2.
! The actual calculation of DeltaNNLO is taken over from the program HDECAY from Michael Spira et al.
! ======================================================================================================
! ============================== Author: Patrick Steppeler =============================================
! ======================================================================================================

!******************
!*** QCD Delta_mb *
!******************

	subroutine DMNLO_CalcDeltaMb(results,gen)

#include "../DMNLO_RenScheme.h"
#include "../DMNLO_Model.h"
#include "../DMNLO_LoopIntegrals.h"

	integer myt, myi, tt, gen
	external DMNLO_CalcInt, DMNLO_ImproveDeltaMb
	double precision DeltaNNLO_QCD, DeltaNNLO_ELW								
	double precision DeltaNLO_At, DeltaNLO_MUE2, DeltaNLO_Ab, DeltaNLO_MUE
	double precision results(5),delMB

	DeltaMb			= 0.0d0			! Total DeltaMb
	DeltaMbQCD		= 0.0d0			! This part gets subtracted in ChiChi2QQ_Couplings.f
	DeltaNNLO_QCD	= 0.0d0			! NNLO QCD Resummation from Spira
	DeltaNNLO_ELW	= 0.0d0			! NNLO ELW Resummation from Spira
	DeltaNLO_Ab		= 0.0d0			! Auxilary Quantity
	DeltaNLO_MUE	= 0.0d0			! " "
	DeltaNLO_At		= 0.0d0			! " "
	DeltaNLO_MUE2	= 0.0d0			! " "

! Extend the Resummation to NNLO or not	
	if((NNLO_Res == 1).AND.(gen == 3)) then
		call DMNLO_ImproveDeltaMb(DeltaNNLO_QCD, DeltaNNLO_ELW)
		Res_Style = 1
	endif

! ------- Gluino-sbottom contribution -------

	call DMNLO_CalcInt(MSf(1,4,gen),MSf(2,4,gen),MGl)
	DeltaNLO_Ab = CF*Alfas/pi*MGl/2d0*(Af(4,gen))*Iabc
	DeltaNLO_MUE = CF*Alfas/pi*MGl/2d0*( - MUE*TB)*Iabc

****************************************************************************************************************************
*** Note Patrick: Here are the different options to include the A_b term. Res_Style = 0 is the old DMNLO variant, while	 ***
*** Res_Style = 1 corresponds to the new Spira variant in agreement with DOI: 10.1103/PhysRevD.68.115001 (See especially ***
*** equations (14), (21) and the text after equation (21)).	Res_Style = 1 is also in agreement with HDECAY.			 ***
***																				 ***
*** Second Note Patrick: DeltaMbQCD is later subtracted in mycouplings.f and contains only the NLO QCD Corrections. It	 ***
*** is important to seperate the NNLO QCD Corrections from this quantity.								 ***
****************************************************************************************************************************

	if(Res_Style == 0) then
		DeltaMbQCD = DeltaMbQCD + DeltaNLO_Ab + DeltaNLO_MUE
		DeltaMb = DeltaMb + DeltaMbQCD
	else
		DeltaMbQCD = DeltaMbQCD + DeltaNLO_MUE/(1 + DeltaNLO_Ab)
		DeltaMb = (DeltaMb  + DeltaNLO_MUE*(1 + Alfas/pi*DeltaNNLO_QCD))/(1 + DeltaNLO_Ab)		
	endif

	results(1) = DeltaNLO_MUE/(1 + DeltaNLO_Ab)
	results(2) = DeltaNLO_MUE*(Alfas/pi*DeltaNNLO_QCD)/(1 + DeltaNLO_Ab)

! ------- Chargino-stop contribution ------- 				
! Takes the real part of mu for CalcInt to avoid problems

	call DMNLO_CalcInt(dreal(MUE),MSf(1,3,gen),MSf(2,3,gen))
	DeltaNLO_At = -(EL*Mf(3,gen)/(sqrt2*MW*SW*SB))**2*1d0/(16*pi**2)*Af(3,gen)*MUE*TB*Iabc
	DeltaNLO_MUE2 = (EL*Mf(3,gen)/(sqrt2*MW*SW*SB))**2*1d0/(16*pi**2)*MUE**2*Iabc

****************************************************************************************************************************
*** Note Patrick: Here are the different options for the electroweak resummation. Res_Style = 0 is again the old DMNLO   ***
*** variant, while Res_Style = 1 corresponds to the new Spira variant in agreement with JHEP06(2011)084 (See especially  ***
*** equations (2.2) and (3.14)). Res_Style = 1 is in agreement with HDECAY. The additional terms for Res_Style = 0  are  ***
*** probably steming from the comparison with micrOMEGAs (See http://arxiv.org/abs/hep-ph/0405253v1 equation (B.66)).    ***
****************************************************************************************************************************

	if(Res_Style ==0) then
		DeltaMb = DeltaMb + DeltaNLO_MUE2 + DeltaNLO_At
		do myt=3,4
			do myi=1,2
				call DMNLO_CalcInt(MSf(myi,myt,gen),M_2,dreal(MUE))
				DeltaMb = DeltaMb + 1d0/(myt-2)*EL**2/SW**2*1d0/(16*pi**2)*MUE*M_2*TB*Usf(1,myi,myt,gen)**2*Iabc
			enddo
		enddo
	else
		DeltaMb = DeltaMb + DeltaNLO_At*(1 + Alfas/pi*DeltaNNLO_ELW)
	endif

	results(3) = DeltaNLO_At
	results(4) = DeltaNLO_At*Alfas/pi*DeltaNNLO_ELW
	results(5) = DeltaMb


	end

!*******************
!*** Integral Iabc *
!*******************


	subroutine DMNLO_CalcInt(mxa,mxb,mxc)

#include "../DMNLO_RenScheme.h"

	double precision mxa,mxb,mxc

	Iabc = 1d0/((mxa**2-mxb**2)*(mxb**2-mxc**2)*(mxc**2-mxa**2))*(mxa**2*mxb**2*dlog(mxa**2/mxb**2)
     -  + mxb**2*mxc**2*dlog(mxb**2/mxc**2) + mxc**2*mxa**2*dlog(mxc**2/mxa**2))

	end

!************************
!*** NNLO contributions *
!************************

	subroutine DMNLO_ImproveDeltaMb(DeltaNNLO_QCD, DeltaNNLO_ELW)
	! Interface between DM@NLO and HDECAY for the calculation of DeltaNNLO_QCD and DeltaNNLO_ELW
	
	
#include "../DMNLO_Model.h"
	double precision DeltaNNLO_QCD, DeltaNNLO_ELW, scale, xxt
	double precision fqcd_hdec, felw_hdec
	
c	write(*,*)'Subroutine ImproveDeltaMb wird aufgerufen.'
	
	! Calculation of DeltaNNLO_QCD	
	scale = QScale											! DMNLO Renormalization Scale
c	scale = (MSf(1,4,3) + MSf(2,4,3) + MGl)/3.0d0			! HDECAY Renormatization Scale
	xxt = Af(3,1) - MUE/TB	
	DeltaNNLO_QCD = fqcd_hdec(scale, MTos, MGl, MSf(1,4,3), MSf(2,4,3), MSf(1,3,3), MSf(2,3,3), 
     & MSf(1,3,1), MSf(2,3,1), MSf(1,4,1), MSf(2,4,1), xxt)
	 
	! Calculation of DeltaNNLO_ELW
c	scale = (MSf(1,3,3) + MSf(2,3,3) + dreal(MUE))/3.0d0	! HDECAY Renormalization Scale
	DeltaNNLO_ELW = felw_hdec(scale,dreal(MUE),MGl,MSf(1,4,3), MSf(2,4,3), MSf(1,3,3), MSf(2,3,3),
     & MTos)
	  
	end                    
																								
! ======================================================================================================
! ================== Begin of copied and changed functions from HDECAY =================================
! ======================================================================================================


	double precision function t134p_hdec(am1intern,am2intern,mu)								! From dmb.f l.2 
	! Generic two loop integral
	
	implicit double precision (a-h,o-z)
	double precision m1intern,m2intern,mu,ll1,ll2
	complex*16 cspen
	
	sp(a) = dreal(cspen(dcmplx(a,0.d0)))
	pi = 4*datan(1.d0)
	zeta2 = pi**2/6
	
	if(am1intern.lt.am2intern)then
       m1intern = am1intern
       m2intern = am2intern
	else
       m1intern = am2intern
       m2intern = am1intern
	endif
	ll1 = dlog(mu**2/m1intern**2)
	ll2 = dlog(mu**2/m2intern**2)
	if(m1intern.eq.m2intern)then
       dummy = 7*(m1intern**2+m2intern**2)/2
     .       + m1intern**2*(ll1**2+3*ll1) + m2intern**2*(ll2**2+3*ll2)
     .       - m1intern**2/2*dlog(m1intern**2/m2intern**2)**2
	else
       dummy = 7*(m1intern**2+m2intern**2)/2
     .       + m1intern**2*(ll1**2+3*ll1) + m2intern**2*(ll2**2+3*ll2)
     .       + (m1intern**2-m2intern**2)*(dlog(m1intern**2/m2intern**2)*dlog(1-m1intern**2/m2intern**2)
     .                       + sp(m1intern**2/m2intern**2)-zeta2)
     .       - m1intern**2/2*dlog(m1intern**2/m2intern**2)**2
	endif
	t134p_hdec = dummy/mu**2
	return
	
	end
	
! ------------------------------------------------------------------------------------------------------
	
	
	double precision function t134_hdec(am1intern,am2intern,am3intern,mu)								! From dmb.f l.35
	! Generic two loop integral
	! Causes Problems for am1intern = am2intern + am3intern or am2intern = am1intern + am3intern or am3intern = am1intern + am2intern
	
	implicit double precision (a-h,o-z)
	double precision m1intern,m2intern,m3intern,mu,ll1,ll2,ll3
	complex*16 sp,cspen,lam,ll,cx,cy,phi
	
	sp(cx) = cspen(cx)
	lam(cx,cy) = (1-cx-cy)**2-4*cx*cy
	phi(cx,cy,ll) = (2*cdlog((1+cx-cy-ll)/2)*cdlog((1-cx+cy-ll)/2)
     .              -cdlog(cx)*cdlog(cy)+2*zeta2
     .              -2*sp((1+cx-cy-ll)/2)-2*sp((1-cx+cy-ll)/2))/ll
	eps = 1.d-15
	rim = dcmplx(1.d0,eps)
	pi = 4*datan(1.d0)
	zeta2 = pi**2/6
	m1intern = dmin1(am1intern,am2intern,am3intern)
	m3intern = dmax1(am1intern,am2intern,am3intern)
	
	if(m1intern.eq.am2intern.and.m3intern.eq.am3intern.or.m1intern.eq.am3intern.and.m3intern.eq.am2intern) m2intern = am1intern
	if(m1intern.eq.am1intern.and.m3intern.eq.am3intern.or.m1intern.eq.am3intern.and.m3intern.eq.am1intern) m2intern = am2intern
	if(m1intern.eq.am1intern.and.m3intern.eq.am2intern.or.m1intern.eq.am2intern.and.m3intern.eq.am1intern) m2intern = am3intern
	cx = m1intern**2/m3intern**2*rim
	cy = m2intern**2/m3intern**2*rim
	ll1 = dlog(mu**2/m1intern**2)
	ll2 = dlog(mu**2/m2intern**2)
	ll3 = dlog(mu**2/m3intern**2)
	ll = cdsqrt(lam(cx,cy))
	dummy = 7*(m1intern**2+m2intern**2+m3intern**2)/2
     .      + m1intern**2*(ll1**2+3*ll1) + m2intern**2*(ll2**2+3*ll2)
     .      + m3intern**2*(ll3**2+3*ll3)
     .      +  (m1intern**2-m2intern**2-m3intern**2)/4*dlog(m2intern**2/m3intern**2)**2
     .      + (-m1intern**2+m2intern**2-m3intern**2)/4*dlog(m1intern**2/m3intern**2)**2
     .      + (-m1intern**2-m2intern**2+m3intern**2)/4*dlog(m1intern**2/m2intern**2)**2
     .      + m3intern**2/2*lam(cx,cy)*phi(cx,cy,ll)
	t134_hdec = dummy/mu**2
	return
	
	end
	  
! ------------------------------------------------------------------------------------------------------	


      double precision function fqcd_hdec(scale,amt,amg,amsb1,amsb2,			! From dmb.f l.464
     .                         amst1,amst2,amsu1,amsu2,amsd1,amsd2,xxt)
	 ! Actual Calculation of DeltaNNLO (QCD Part)
	 ! List of Arguments:
	 ! amt -- Top pole mass
	 ! amg, amsb1, amsb2, amst1, amst2, amsu1, amsu2, amsd1, amsd2 -- Loop corrected pole masses of
	 ! Gluino, Sbottom1/2, Stop1/2, Sup1/2, Sdown1/2
	 ! For scale and xxt see main subroutine. The trilinear coupling A_u is in the DRbar or MSbar
	 ! Scheme (Equal for NLO)
	 
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amt,amg,amsb1,amsb2,amst1,amst2,
     .                 amsu1,amsu2,amsd1,amsd2
      double precision mt,mg,mb1,mb2,mt1,mt2,ms1,ms2,mu
      complex*16 sp,cspen,xt,xgl,xt1,xt2,xb1,xb2,xs1,xs2,xq
      double precision m,mq
      double precision anomalous
      double precision a
      complex*16 xp,xm
      sp(r) = cspen(r)
      xp(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . +cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      xm(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . -cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      ca = 3
      cf = 4/3.d0
      tr = 1/2.d0
      nu = 2
      nd = 2
      nf = nu+nd+1

      fnorm = 4/amg**2

      rim = dcmplx(1.d0,eps)

      mq  = amt
      mt  = amt
      mg  = amg
      mu  = amg
      mt1 = amsb1
      mt2 = amsb2
      mb1 = amsb1
      mb2 = amsb2

      xq  = amt**2/amg**2   * rim
      xt  = amt**2/amg**2   * rim
      xgl = amg**2/amg**2   * rim
      xt1 = amsb1**2/amg**2 * rim
      xt2 = amsb2**2/amg**2 * rim
      xb1 = amsb1**2/amg**2 * rim
      xb2 = amsb2**2/amg**2 * rim

      r12=(log(xb1)**2*xb1*(-2*xb1*xb2**2+4*xb1*xb2-2*xb1+xb2**2-2*
     . xb2+1)+4*log(xb1)*xb1*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1
     . )+log(xb2)**2*xb2*(2*xb1**2*xb2-xb1**2-4*xb1*xb2+2*xb1+2*xb2-1
     . )+4*log(xb2)*xb2*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1)+4*
     . t134p(mb1,mb1,mg)*xb1*(xb2**2-2*xb2+1)+4*t134p(mb1,mg,mg)*(-
     . xb1*xb2**2+2*xb1*xb2-xb1-xb2**2+2*xb2-1)+4*t134p(mb2,mb2,mg)*
     . xb2*(-xb1**2+2*xb1-1)+4*t134p(mb2,mg,mg)*(xb1**2*xb2+xb1**2-2*
     . xb1*xb2-2*xb1+xb2+1)+4*t134p(mg,mg,mg)*(-xb1**2+2*xb1+xb2**2-2
     . *xb2))/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*
     . xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2
     . **2-xb2))

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r32=(log(xb1)**2*xb1*(-xb2**2+2*xb2-1)+2*log(xb1)*xb1*(2*xb1*
     . xb2**2-4*xb1*xb2+2*xb1+xb2**2-2*xb2+1)+log(xb2)**2*xb2*(xb1**2
     . -2*xb1+1)+2*log(xb2)*xb2*(-2*xb1**2*xb2-xb1**2+4*xb1*xb2+2*xb1
     . -2*xb2-1)+2*t134p(mb1,mg,mg)*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2
     . **2+2*xb2-1)+2*t134p(mb2,mg,mg)*(-2*xb1**2*xb2+xb1**2+4*xb1*
     . xb2-2*xb1-2*xb2+1)+3*t134p(mg,mg,mg)*(xb1**2*xb2-xb1**2-xb1*
     . xb2**2+xb1+xb2**2-xb2)+14*(xb1**2*xb2-xb1*xb2**2-xb1+xb2))/(8*
     . (xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-
     . 2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r52=(log(xb1)**2*xb1*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2**2+2*xb2
     . -1)+4*log(xb1)*xb1*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+
     . log(xb2)**2*xb2*(-2*xb1**2*xb2+xb1**2+4*xb1*xb2-2*xb1-2*xb2+1)
     . +4*log(xb2)*xb2*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+2*
     . t134p(mb1,mb2,mg)*(xb1**2*xb2-xb1**2-xb1*xb2**2+xb1+xb2**2-xb2
     . )+2*t134p(mb1,mg,mg)*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1
     . )+2*t134p(mb2,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1)
     . )/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**
     . 2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-
     . xb2))

      r62=(log(xb1)**2*xb1*(xb2**2-2*xb2+1)+2*log(xb1)*xb1*(-xb2**2+2
     . *xb2-1)+log(xb2)**2*xb2*(-xb1**2+2*xb1-1)+2*log(xb2)*xb2*(xb1
     . **2-2*xb1+1)+2*t134p(mb1,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1+xb2
     . **2-2*xb2+1)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2-xb1**2+2*xb1*xb2+
     . 2*xb1-xb2-1)+2*t134p(mg,mg,mg)*(xb1**2*xb2+xb1**2-xb1*xb2**2-3
     . *xb1-xb2**2+3*xb2))/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1
     . **2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1
     . -xb2**3+2*xb2**2-xb2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      r82=(log(xb1)**2*(xb1*xb2**2-2*xb1*xb2+xb1-2*xb2**2+4*xb2-2)+4*
     . log(xb1)*(-xb1**2*xb2**2+2*xb1**2*xb2-xb1**2+xb1*xb2**2-2*xb1*
     . xb2+xb1+xb2**2-2*xb2+1)+log(xb2)**2*(-xb1**2*xb2+2*xb1**2+2*
     . xb1*xb2-4*xb1-xb2+2)+4*log(xb2)*(xb1**2*xb2**2-xb1**2*xb2-xb1
     . **2-2*xb1*xb2**2+2*xb1*xb2+2*xb1+xb2**2-xb2-1)+2*t134p(mb1,mg,
     . mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+2*t134p(mb2,mg,
     . mg)*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+10*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2))/(8*(xb1**3*xb2**2-2*xb1**3*
     . xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*
     . xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      mt1 = amst1
      mt2 = amst2
      xt1 = amst1**2/amg**2 * rim
      xt2 = amst2**2/amg**2 * rim

      ralsca2=(-((9*log(xb1)-10)*(xb2-1)*log(xb1)*xb1-(9*log(xb2)-10)
     . *(xb1-1)*log(xb2)*xb2))/(48*(xb1-xb2)*(xb1-1)*(xb2-1))

      ralscf2=(-((xb1-1)*log(xb2)*xb2-(xb2-1)*log(xb1)*xb1))/(8*(xb1-
     . xb2)*(xb1-1)*(xb2-1))

      rmb12=(-((((2*xb1**3-3*xb2-(xb2-6)*xb1**2-2*(xb2+1)*xb1)*xb1+(
     . xb1**2-xb2)*(xb1-1)**2*log(-(xb1-1)))*(xb2-1)+(xb1-1)**2*log(
     . xb2)*xb1**2*xb2)*log(xb1)-(((xb1**2-xb2)*(xb1+1)*(xb2-1)*log(
     . xb1)**2-(xb1-1)**2*log(xb2)**2*xb2)*xb1+((xb1-xb2)*(xb2-1)+(
     . xb1-1)*log(xb2)*xb2)*((xb1+5)*xb1+(xb1-1)**2*log(-(xb1-1)))*(
     . xb1-1))))/(4*(xb1-xb2)**2*(xb1-1)**2*(xb2-1)*xb1)

      rmb22=(-(((xb1-xb2**2)*(xb1-1)*(xb2+1)*log(xb2)**2+(xb2-1)**2*
     . log(xb1)**2*xb1)*xb2+((xb2+5)*xb2+(xb2-1)**2*log(-(xb2-1)))*(
     . xb1-xb2)*(xb1-1)*(xb2-1)-((xb2-1)**2*log(-(xb2-1))-log(xb2)*
     . xb2**2+(xb2+5)*xb2)*(xb2-1)**2*log(xb1)*xb1+((2*(xb2**2+3*xb2-
     . 1)*xb2-(xb2**2+2*xb2+3)*xb1)*xb2-(xb1-xb2**2)*(xb2-1)**2*log(-
     . (xb2-1)))*(xb1-1)*log(xb2)))/(4*(xb1-xb2)**2*(xb1-1)*(xb2-1)**
     . 2*xb2)

      rmgca2=(-((3*log(xb1)-14)*(xb1+1)*(xb2-1)**2*log(xb1)*xb1-(3*
     . log(xb2)-14)*(xb1-1)**2*(xb2+1)*log(xb2)*xb2-28*(xb1-xb2)*(xb1
     . -1)*(xb2-1)))/(16*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsu1
      ms2 = amsu2
      xs1 = amsu1**2/amg**2 * rim
      xs2 = amsu2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      ru102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rualstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rumgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsd1
      ms2 = amsd2
      xs1 = amsd1**2/amg**2 * rim
      xs2 = amsd2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      rd102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rdalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rdmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsb1
      ms2 = amsb2
      xs1 = amsb1**2/amg**2 * rim
      xs2 = amsb2**2/amg**2 * rim

      ans1=log(xb1)**2*xb1*(xb1**2*xb2-xb1**2-4*xb1*xb2**2+6*xb1*xb2-
     . 2*xb1-xb2**3+3*xb2**2-2*xb2)+2*log(xb1)*log(xb2)*xb1*xb2*(xb1
     . **2-2*xb1-xb2**2+2*xb2)+4*log(xb1)*xb1*(-2*xb1**2*xb2+2*xb1**2
     . +3*xb1*xb2**2-2*xb1*xb2-xb1+xb2**3-4*xb2**2+3*xb2)+log(xb2)**2
     . *xb2*(xb1**3+4*xb1**2*xb2-3*xb1**2-xb1*xb2**2-6*xb1*xb2+2*xb1+
     . xb2**2+2*xb2)+4*log(xb2)*xb2*(-xb1**3-3*xb1**2*xb2+4*xb1**2+2*
     . xb1*xb2**2+2*xb1*xb2-3*xb1-2*xb2**2+xb2)+2*t134p(mb1,mb2,mg)*(
     . -xb1**3+xb1**2*xb2+2*xb1**2-xb1*xb2**2-2*xb1+xb2**3-2*xb2**2+2
     . *xb2)+2*t134p(mb1,mg,mg)*(xb1**3-2*xb1**2*xb2-xb1**2+xb1*xb2**
     . 2+2*xb1*xb2-xb2**2)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2+xb1**2+2*
     . xb1*xb2**2-2*xb1*xb2-xb2**3+xb2**2)+12*(xb1**3*xb2-xb1**3-xb1*
     . xb2**3+xb1+xb2**3-xb2)
      rb102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rbalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rbmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amst1
      ms2 = amst2
      xs1 = amst1**2/amg**2 * rim
      xs2 = amst2**2/amg**2 * rim

      ans14=-8*log(xs2)*xb1*xb2**2*xs2-4*log(xs2)*xb1*xq**2*xs2-4*log
     . (xs2)*xb1*xq*xs2**2-4*log(xs2)*xb1*xq*xs2+8*log(xs2)*xb1*xs2**
     . 3-16*log(xs2)*xb1*xs2**2+8*log(xs2)*xb1*xs2-4*log(xs2)*xb2**2*
     . xq**2*xs2-4*log(xs2)*xb2**2*xq*xs2**2-4*log(xs2)*xb2**2*xq*xs2
     . +8*log(xs2)*xb2**2*xs2**3-16*log(xs2)*xb2**2*xs2**2+8*log(xs2)
     . *xb2**2*xs2+4*log(xs2)*xb2*xq**2*xs2+4*log(xs2)*xb2*xq*xs2**2+
     . 4*log(xs2)*xb2*xq*xs2-8*log(xs2)*xb2*xs2**3+16*log(xs2)*xb2*
     . xs2**2-8*log(xs2)*xb2*xs2
      ans13=-log(xs2)**2*xb1*xb2**2*xq**2*xs2+log(xs2)**2*xb1*xb2**2*
     . xs2**3-2*log(xs2)**2*xb1*xb2**2*xs2**2+log(xs2)**2*xb1*xb2**2*
     . xs2+log(xs2)**2*xb1*xq**2*xs2-log(xs2)**2*xb1*xs2**3+2*log(xs2
     . )**2*xb1*xs2**2-log(xs2)**2*xb1*xs2+log(xs2)**2*xb2**2*xq**2*
     . xs2-log(xs2)**2*xb2**2*xs2**3+2*log(xs2)**2*xb2**2*xs2**2-log(
     . xs2)**2*xb2**2*xs2-log(xs2)**2*xb2*xq**2*xs2+log(xs2)**2*xb2*
     . xs2**3-2*log(xs2)**2*xb2*xs2**2+log(xs2)**2*xb2*xs2-4*log(xs2)
     . *xb1**2*xb2*xq**2*xs2-4*log(xs2)*xb1**2*xb2*xq*xs2**2-4*log(
     . xs2)*xb1**2*xb2*xq*xs2+8*log(xs2)*xb1**2*xb2*xs2**3-16*log(xs2
     . )*xb1**2*xb2*xs2**2+8*log(xs2)*xb1**2*xb2*xs2+4*log(xs2)*xb1**
     . 2*xq**2*xs2+4*log(xs2)*xb1**2*xq*xs2**2+4*log(xs2)*xb1**2*xq*
     . xs2-8*log(xs2)*xb1**2*xs2**3+16*log(xs2)*xb1**2*xs2**2-8*log(
     . xs2)*xb1**2*xs2+4*log(xs2)*xb1*xb2**2*xq**2*xs2+4*log(xs2)*xb1
     . *xb2**2*xq*xs2**2+4*log(xs2)*xb1*xb2**2*xq*xs2-8*log(xs2)*xb1*
     . xb2**2*xs2**3+16*log(xs2)*xb1*xb2**2*xs2**2+ans14
      ans12=-4*log(xb2)**2*xb1*xb2*xq*xs1-8*log(xb2)**2*xb1*xb2*xq*
     . xs2**2+4*log(xb2)**2*xb1*xb2*xq*xs2-4*log(xb2)**2*xb1*xb2*xq+2
     . *log(xb2)**2*xb1*xb2*xs1*xs2**2-4*log(xb2)**2*xb1*xb2*xs1*xs2+
     . 2*log(xb2)**2*xb1*xb2*xs1+2*log(xb2)**2*xb1*xb2*xs2**3-4*log(
     . xb2)**2*xb1*xb2*xs2**2+2*log(xb2)**2*xb1*xb2*xs2+2*log(xb2)**2
     . *xb2*xq**3-log(xb2)**2*xb2*xq**2*xs1-5*log(xb2)**2*xb2*xq**2*
     . xs2-4*log(xb2)**2*xb2*xq**2+2*log(xb2)**2*xb2*xq*xs1*xs2+2*log
     . (xb2)**2*xb2*xq*xs1+4*log(xb2)**2*xb2*xq*xs2**2-2*log(xb2)**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb2*xq-log(xb2)**2*xb2*xs1*xs2**2+2*
     . log(xb2)**2*xb2*xs1*xs2-log(xb2)**2*xb2*xs1-log(xb2)**2*xb2*
     . xs2**3+2*log(xb2)**2*xb2*xs2**2-log(xb2)**2*xb2*xs2+log(xs2)**
     . 2*xb1**2*xb2*xq**2*xs2-log(xs2)**2*xb1**2*xb2*xs2**3+2*log(xs2
     . )**2*xb1**2*xb2*xs2**2-log(xs2)**2*xb1**2*xb2*xs2-log(xs2)**2*
     . xb1**2*xq**2*xs2+log(xs2)**2*xb1**2*xs2**3-2*log(xs2)**2*xb1**
     . 2*xs2**2+log(xs2)**2*xb1**2*xs2+ans13
      ans11=4*log(xb1)**2*xb1*xq**2-2*log(xb1)**2*xb1*xq*xs1*xs2-2*
     . log(xb1)**2*xb1*xq*xs1-4*log(xb1)**2*xb1*xq*xs2**2+2*log(xb1)
     . **2*xb1*xq*xs2-2*log(xb1)**2*xb1*xq+log(xb1)**2*xb1*xs1*xs2**2
     . -2*log(xb1)**2*xb1*xs1*xs2+log(xb1)**2*xb1*xs1+log(xb1)**2*xb1
     . *xs2**3-2*log(xb1)**2*xb1*xs2**2+log(xb1)**2*xb1*xs2+2*log(xb2
     . )**2*xb1**2*xb2*xq**3-log(xb2)**2*xb1**2*xb2*xq**2*xs1-5*log(
     . xb2)**2*xb1**2*xb2*xq**2*xs2-4*log(xb2)**2*xb1**2*xb2*xq**2+2*
     . log(xb2)**2*xb1**2*xb2*xq*xs1*xs2+2*log(xb2)**2*xb1**2*xb2*xq*
     . xs1+4*log(xb2)**2*xb1**2*xb2*xq*xs2**2-2*log(xb2)**2*xb1**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb1**2*xb2*xq-log(xb2)**2*xb1**2*xb2*
     . xs1*xs2**2+2*log(xb2)**2*xb1**2*xb2*xs1*xs2-log(xb2)**2*xb1**2
     . *xb2*xs1-log(xb2)**2*xb1**2*xb2*xs2**3+2*log(xb2)**2*xb1**2*
     . xb2*xs2**2-log(xb2)**2*xb1**2*xb2*xs2-4*log(xb2)**2*xb1*xb2*xq
     . **3+2*log(xb2)**2*xb1*xb2*xq**2*xs1+10*log(xb2)**2*xb1*xb2*xq
     . **2*xs2+8*log(xb2)**2*xb1*xb2*xq**2-4*log(xb2)**2*xb1*xb2*xq*
     . xs1*xs2+ans12
      ans10=5*log(xb1)**2*xb1*xb2**2*xq**2*xs2+4*log(xb1)**2*xb1*xb2
     . **2*xq**2-2*log(xb1)**2*xb1*xb2**2*xq*xs1*xs2-2*log(xb1)**2*
     . xb1*xb2**2*xq*xs1-4*log(xb1)**2*xb1*xb2**2*xq*xs2**2+2*log(xb1
     . )**2*xb1*xb2**2*xq*xs2-2*log(xb1)**2*xb1*xb2**2*xq+log(xb1)**2
     . *xb1*xb2**2*xs1*xs2**2-2*log(xb1)**2*xb1*xb2**2*xs1*xs2+log(
     . xb1)**2*xb1*xb2**2*xs1+log(xb1)**2*xb1*xb2**2*xs2**3-2*log(xb1
     . )**2*xb1*xb2**2*xs2**2+log(xb1)**2*xb1*xb2**2*xs2+4*log(xb1)**
     . 2*xb1*xb2*xq**3-2*log(xb1)**2*xb1*xb2*xq**2*xs1-10*log(xb1)**2
     . *xb1*xb2*xq**2*xs2-8*log(xb1)**2*xb1*xb2*xq**2+4*log(xb1)**2*
     . xb1*xb2*xq*xs1*xs2+4*log(xb1)**2*xb1*xb2*xq*xs1+8*log(xb1)**2*
     . xb1*xb2*xq*xs2**2-4*log(xb1)**2*xb1*xb2*xq*xs2+4*log(xb1)**2*
     . xb1*xb2*xq-2*log(xb1)**2*xb1*xb2*xs1*xs2**2+4*log(xb1)**2*xb1*
     . xb2*xs1*xs2-2*log(xb1)**2*xb1*xb2*xs1-2*log(xb1)**2*xb1*xb2*
     . xs2**3+4*log(xb1)**2*xb1*xb2*xs2**2-2*log(xb1)**2*xb1*xb2*xs2-
     . 2*log(xb1)**2*xb1*xq**3+log(xb1)**2*xb1*xq**2*xs1+5*log(xb1)**
     . 2*xb1*xq**2*xs2+ans11
      ans9=-4*t134(mb1,mq,ms2,mg)*xb1*xb2+2*t134(mb1,mq,ms2,mg)*xb1*
     . xq**2-4*t134(mb1,mq,ms2,mg)*xb1*xq*xs2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xq+2*t134(mb1,mq,ms2,mg)*xb1*xs2**2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xs2+2*t134(mb1,mq,ms2,mg)*xb1+2*t134(mb1,mq,ms2,mg)*xb2**2
     . *xq**3-6*t134(mb1,mq,ms2,mg)*xb2**2*xq**2*xs2-4*t134(mb1,mq,
     . ms2,mg)*xb2**2*xq**2+6*t134(mb1,mq,ms2,mg)*xb2**2*xq*xs2**2+2*
     . t134(mb1,mq,ms2,mg)*xb2**2*xq-2*t134(mb1,mq,ms2,mg)*xb2**2*xs2
     . **3+4*t134(mb1,mq,ms2,mg)*xb2**2*xs2**2-2*t134(mb1,mq,ms2,mg)*
     . xb2**2*xs2-4*t134(mb1,mq,ms2,mg)*xb2*xq**3+12*t134(mb1,mq,ms2,
     . mg)*xb2*xq**2*xs2+8*t134(mb1,mq,ms2,mg)*xb2*xq**2-12*t134(mb1,
     . mq,ms2,mg)*xb2*xq*xs2**2-4*t134(mb1,mq,ms2,mg)*xb2*xq+4*t134(
     . mb1,mq,ms2,mg)*xb2*xs2**3-8*t134(mb1,mq,ms2,mg)*xb2*xs2**2+4*
     . t134(mb1,mq,ms2,mg)*xb2*xs2+2*t134(mb1,mq,ms2,mg)*xq**3-6*t134
     . (mb1,mq,ms2,mg)*xq**2*xs2-4*t134(mb1,mq,ms2,mg)*xq**2+6*t134(
     . mb1,mq,ms2,mg)*xq*xs2**2+2*t134(mb1,mq,ms2,mg)*xq-2*t134(mb1,
     . mq,ms2,mg)*xs2**3+4*t134(mb1,mq,ms2,mg)*xs2**2-2*t134(mb1,mq,
     . ms2,mg)*xs2-2*log(xb1)**2*xb1*xb2**2*xq**3+log(xb1)**2*xb1*xb2
     . **2*xq**2*xs1+ans10
      ans8=4*t134(mb1,mq,ms1,mg)*xb2*xq**2*xs1+8*t134(mb1,mq,ms1,mg)*
     . xb2*xq**2*xs2+8*t134(mb1,mq,ms1,mg)*xb2*xq**2-8*t134(mb1,mq,
     . ms1,mg)*xb2*xq*xs1*xs2-8*t134(mb1,mq,ms1,mg)*xb2*xq*xs1-4*t134
     . (mb1,mq,ms1,mg)*xb2*xq*xs2**2+8*t134(mb1,mq,ms1,mg)*xb2*xq*xs2
     . -4*t134(mb1,mq,ms1,mg)*xb2*xq+4*t134(mb1,mq,ms1,mg)*xb2*xs1*
     . xs2**2-8*t134(mb1,mq,ms1,mg)*xb2*xs1*xs2+4*t134(mb1,mq,ms1,mg)
     . *xb2*xs1+2*t134(mb1,mq,ms1,mg)*xq**3-2*t134(mb1,mq,ms1,mg)*xq
     . **2*xs1-4*t134(mb1,mq,ms1,mg)*xq**2*xs2-4*t134(mb1,mq,ms1,mg)*
     . xq**2+4*t134(mb1,mq,ms1,mg)*xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*
     . xq*xs1+2*t134(mb1,mq,ms1,mg)*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xq*xs2+2*t134(mb1,mq,ms1,mg)*xq-2*t134(mb1,mq,ms1,mg)*xs1*xs2
     . **2+4*t134(mb1,mq,ms1,mg)*xs1*xs2-2*t134(mb1,mq,ms1,mg)*xs1+2*
     . t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2**2*xq*xs2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq+2*t134(mb1,
     . mq,ms2,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*
     . xs2+2*t134(mb1,mq,ms2,mg)*xb1*xb2**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2*xq**2+8*t134(mb1,mq,ms2,mg)*xb1*xb2*xq*xs2+8*t134(mb1,mq,
     . ms2,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms2,mg)*xb1*xb2*xs2**2+8*t134
     . (mb1,mq,ms2,mg)*xb1*xb2*xs2+ans9
      ans7=2*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms1,
     . mg)*xb1*xb2**2*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq+2*
     . t134(mb1,mq,ms1,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xb1*xb2**2*xs2+2*t134(mb1,mq,ms1,mg)*xb1*xb2**2-4*t134(mb1,mq,
     . ms1,mg)*xb1*xb2*xq**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xq*xs2+8*
     . t134(mb1,mq,ms1,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms1,mg)*xb1*xb2*
     . xs2**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xs2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xb2+2*t134(mb1,mq,ms1,mg)*xb1*xq**2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xq+2*t134(mb1,mq,ms1,mg)
     . *xb1*xs2**2-4*t134(mb1,mq,ms1,mg)*xb1*xs2+2*t134(mb1,mq,ms1,mg
     . )*xb1+2*t134(mb1,mq,ms1,mg)*xb2**2*xq**3-2*t134(mb1,mq,ms1,mg)
     . *xb2**2*xq**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2**2*xq**2*xs2-4*
     . t134(mb1,mq,ms1,mg)*xb2**2*xq**2+4*t134(mb1,mq,ms1,mg)*xb2**2*
     . xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs1+2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs2+2
     . *t134(mb1,mq,ms1,mg)*xb2**2*xq-2*t134(mb1,mq,ms1,mg)*xb2**2*
     . xs1*xs2**2+4*t134(mb1,mq,ms1,mg)*xb2**2*xs1*xs2-2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2*xq**3+ans8
      ans15=(xq**2-2*xq*xs1-2*xq+xs1**2-2*xs1+1)
      ans6=ans7*ans15
      ans5=-ans6
      ans4=(4*(2*(xs1-1)**2-xq**2-(xs1+1)*xq)-(xs1-1+xq)*(xs1-1-xq)*
     . log(xs1))*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*(xb1-1)*(
     . xb2-1)*log(xs1)*xs1+ans5
      ans3=-ans4
      ans17=4*(3*((xs1+1)*xs1+xs2**2+xs2-1)+2*(3*xq**2+2)*xq+((2*xq**
     . 2-5*xq+8)*xq-3*(xs1+xs2)**2)*(xs1+xs2)-(4*xq-3*xs1*xs2)*(2*xs1
     . **2+3*xs1*xs2+2*xs2**2)+2*(3*xs1**2-8*xs1*xs2+3*xs2**2)*(xs1+
     . xs2)*xq-(13*xq**2+4*xs1*xs2)*xq**2+((3*xs1**2-2*xs1*xs2+3*xs2
     . **2)*xq-2*(xs1+xs2)*xq**2-3*(xs1+xs2)*xs1*xs2+6*xq**3)*(xq-xs1
     . )*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)+4*(2*((xs2-1)**2+xs1**2+
     . 2*(2*xs2-1)*xs1+(5*xq**2-1)*xq)-(xs1**2-4*xs1*xs2+xs2**2)*xq-2
     . *(xq**2+2*xs1*xs2)*(xs1+xs2)+(xs1-1-xq)*(2*(xs2+1)*xq-(xs2-1)
     . **2-xq**2)*log(xs1)*xs1+(2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xs2-1
     . -xq)*log(xs2)*xs2+(xs1+xs2-6*xq)*xq+((xs1+xs2)*xq-4*xq**2+2*
     . xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq)*
     . xq
      ans16=-2*(2*((xs1-1)**2*xs1-xq**3-(3*xs1-1)*xq*xs1+(3*xs1+1)*xq
     . **2)+(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1)*xb2-(2*((xs1+1
     . )*xq-(xs1-1)**2)*xb2-(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1
     . ))*xb1)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*t134(mg,mq,
     . ms1,mg)-2*((xs2-1)**2+xs1**2+2*(2*xs2-1)*xs1+2*(xq+1)*(xq-1)*
     . xq+(xq-2*xs1*xs2)*(xs1+xs2)-(xs1**2-4*xs1*xs2+xs2**2)*xq-(xq**
     . 2-xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq
     . )**2*xq-2*((2*((xs2-1)**2*xs2-xq**3-(3*xs2-1)*xq*xs2+(3*xs2+1)
     . *xq**2)+(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-xs2-1)*xb2-(2*((
     . xs2+1)*xq-(xs2-1)**2)*xb2-(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-
     . xs2-1))*xb1)*(xb1-xb2)*t134(mg,mq,ms2,mg)+((xq-xs1+xb2)*t134(
     . mb2,mq,ms1,mg)+(xq-xs2+xb2)*t134(mb2,mq,ms2,mg))*(2*(xs2+1)*xq
     . -(xs2-1)**2-xq**2)*(xb1-1)**2)*(2*(xs1+1)*xq-(xs1-1)**2-xq**2)
     . +ans17
      ans2=2*(2*(xs1+xs2-2*xq)-log(xs2)*xs2-log(xs1)*xs1+2*log(xq)*xq
     . )*((xb1-1)**2*log(xb2)*xb2-(xb2-1)**2*log(xb1)*xb1)*(2*(xs1+1)
     . *xq-(xs1-1)**2-xq**2)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)+ans3+
     . ans16
      ans1=-ans2
      rt102=ans1/(8*((xs1-1)**2+xq**2-2*(xs1+1)*xq)*((xs2-1)**2+xq**2-
     . 2*(xs2+1)*xq)*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      rtalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rtmgtr2=((xs2-6+xs1-2*xq+(xq+1)*log(xs2)+(xq+1)*log(xs1)-2*log(
     . xq)*xq-(xs2-1-xq)*log(xp(ms2))*xp(ms2)-(xs1-1-xq)*log(xp(ms1))
     . *xp(ms1)-(xs2-1-xq)*log(xm(ms2))*xm(ms2)-(xs1-1-xq)*log(xm(ms1
     . ))*xm(ms1)+(xs2-1-xq)*log(xp(ms2)-1)*xp(ms2)+(xs1-1-xq)*log(xp
     . (ms1)-1)*xp(ms1)+(xs2-1-xq)*log(xm(ms2)-1)*xm(ms2)+(xs1-1-xq)*
     . log(xm(ms1)-1)*xm(ms1))*((xb1+1)*(xb2-1)**2*log(xb1)*xb1-(xb1-
     . 1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs1)+log(xs2))*(xq+1)-2*
     . log(xq)*xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*((log(xm(ms2)
     . -1)-log(xm(ms2)))*xm(ms2)+(log(xp(ms2)-1)-log(xp(ms2)))*xp(ms2
     . ))*(xs2-1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*((log(xm(ms1)-1)-log
     . (xm(ms1)))*xm(ms1)+(log(xp(ms1)-1)-log(xp(ms1)))*xp(ms1))*(xs1
     . -1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*(xs2-6+xs1-2*xq)*(xb1-xb2)*
     . (xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ans2=(((3*xs1-1)*(xs1-1)+3*xq**2-2*(3*xs1+2)*xq-(xq-xs1+1)*(xq-
     . xs1-1)*xb2-(((xs1-1)*(xs1-3)+xq**2-2*(xs1+2)*xq)*xb2+(xq-xs1+1
     . )*(xq-xs1-1))*xb1)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*t134(mg,mq,
     . ms1,mg)-2*(2*(2*xs2-1+2*xq-(xq-xs2)**2)*log(xs1)*xs1-3*((xs2+2
     . +xs1)*xq+(xs1-1)*(xs2-1)-3*xq**2)*(xs1-xs2)-2*(2*xs1-1+2*xq-(
     . xq-xs1)**2)*log(xs2)*xs2)*(xb1-1)*(xb2-1)*log(xq)*xq)*(xb1-xb2
     . )
      ans1=((2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xs2-1+xq)*(log(xs2)-6)*
     . log(xs2)*xs2-(xs1-1+xq)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(log(
     . xs1)-6)*log(xs1)*xs1)*(xb1-xb2)*(xb1-1)*(xb2-1)+((t134(mb1,mq,
     . ms1,mg)-t134(mb1,mq,ms2,mg))*(xb1+1)*(xb2-1)**2-(t134(mb2,mq,
     . ms1,mg)-t134(mb2,mq,ms2,mg))*(xb1-1)**2*(xb2+1))*(2*(xs1+1)*xq
     . -(xs1-1)**2-xq**2)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)-((xs2+2+xs1
     . )*xq+(xs1-1)*(xs2-1)-3*xq**2)*(xb1-xb2)*(xb1-1)*(xb2-1)*(xs1-
     . xs2)*log(xq)**2*xq+14*(xs2-1+xs1-(xq+xs1)*(xq+xs2)+2*(xq**2-
     . xs1*xs2)*xq)*(xb1-xb2)*(xb1-1)*(xb2-1)*(xs1-xs2)-((3*xs2-1)*(
     . xs2-1)+3*xq**2-2*(3*xs2+2)*xq-(xq-xs2+1)*(xq-xs2-1)*xb2-(((xs2
     . -1)*(xs2-3)+xq**2-2*(xs2+2)*xq)*xb2+(xq-xs2+1)*(xq-xs2-1))*xb1
     . )*(2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xb1-xb2)*t134(mg,mq,ms2,mg)
     . +ans2
      rtp102=ans1/(2*((xs1-1)**2+xq**2-2*(xs1+1)*xq)*((xs2-1)**2+xq**2-
     . 2*(xs2+1)*xq)*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2*(xs1-xs2))

      rmgtrt2=(-((log(xp(ms2)-1)-log(xp(ms2)))*xp(ms2)+log(xs1)-log(
     . xs2)-(log(xp(ms1)-1)-log(xp(ms1)))*xp(ms1)+(log(xm(ms2)-1)-log
     . (xm(ms2)))*xm(ms2)-(log(xm(ms1)-1)-log(xm(ms1)))*xm(ms1))*((
     . xb1+1)*(xb2-1)**2*log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2
     . +2*(xb1-xb2)*(xb1-1)*(xb2-1)))/(4*(xb1-xb2)*(xb1-1)**2*(xb2-1)
     . **2)

      fact = amt**2*xxt/amg**3
      rtp102 = fact * rtp102
      fact = 2*amt**2*xxt/(amst1**2-amst2**2)/amg
      rmgtrt2 = fact * rmgtrt2
      rtp102 = rtp102 + rmgtrt2

      rctca = ralsca2 + rmgca2
      rctcf = ralscf2 + rmb12 + rmb22
      ructtr = rualstr2 + rumgtr2
      rdcttr = rdalstr2 + rdmgtr2
      rbcttr = rbalstr2 + rbmgtr2
      rtcttr = rtalstr2 + rtmgtr2

      bo  =-2*fi(amsb1**2,amsb2**2,amg**2)

      rca = r12 + r22/4 + r32 + r52 + r62 + rctca
      rcf = -r22/2 - r42/2 - 2*r52 - r72/2 + r82 - r92/2 + rctcf
      rtr = nu*(ru102 + ructtr) + nd*(rd102 + rdcttr)
      rtrb= rb102 + rbcttr
      rtrt= rt102 + rtp102 + rtcttr + bo/3*dlog(amg**2/amt**2)/fnorm

      rqcd = ca*rca + cf*rcf + tr*(rtr+rtrb+rtrt)
      rqcd = rqcd/bo
      anomalous = - cf/4
      finscale = (11*ca-4*tr*nf)/12*dlog(scale**2/amg**2)
      fqcd_hdec = dreal(rqcd)*fnorm + anomalous + finscale

      return
	  
      end
	  
! ------------------------------------------------------------------------------------------------------

	double precision function felw_hdec(scale,amu,amg,amsb1,amsb2,amst1,amst2,amt)	! From dmb.f l.464
	 ! Actual Calculation of DeltaNNLO (ELW Part)
	 ! List of Arguments:
	 ! amt -- Top pole mass
	 ! amu -- Higgs Parameter mue
	 ! amg, amsb1, amsb2, amst1, amst2 -- Loop corrected pole masses of Gluino, Sbottom1/2, Stop1/2
	 ! For scale see main subroutine.
	 
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amg,amsb1,amsb2,amst1,amst2,amu,amt
      double precision mg,mb1,mb2,mt1,mt2,mu,mt
      double precision anomalous
      double precision a
      complex*16 sp,cpen,xgl,xt1,xt2,xt,xb1,xb2
      sp(r) = cspen(r)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      fnorm = 1/amu**2

      cf = 4/3.d0

      rim = dcmplx(1.d0,eps)

      mt  = amt
      mu  = dabs(amu)
      mg  = amg
      mt1 = amst1
      mt2 = amst2
      mb1 = amsb1
      mb2 = amsb2

      xt  = amt**2/amu**2   * rim
      xgl = amg**2/amu**2   * rim
      xt1 = amst1**2/amu**2 * rim
      xt2 = amst2**2/amu**2 * rim
      xb1 = amsb1**2/amu**2 * rim
      xb2 = amsb2**2/amu**2 * rim

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      ans3=t134(mt2,mg,mt,mg)*xgl*(-xb1*xgl*xt1+xb1*xgl-xb1*xt*xt1+
     . xb1*xt+xb1*xt1*xt2-xb1*xt2-xb2*xgl*xt1+xb2*xgl-xb2*xt*xt1+xb2*
     . xt+xb2*xt1*xt2-xb2*xt2+2*xgl**2*xt1-2*xgl**2+2*xgl*xt*xt1-2*
     . xgl*xt-2*xgl*xt1*xt2+2*xgl*xt2)+t134(mu,mb1,mt,mg)*xgl*(-xb1*
     . xb2*xt1+xb1*xb2*xt2+xb1*xgl*xt1-xb1*xgl*xt2-xb2*xt*xt1+xb2*xt*
     . xt2+xb2*xt1-xb2*xt2+xgl*xt*xt1-xgl*xt*xt2-xgl*xt1+xgl*xt2)+
     . t134(mu,mb2,mt,mg)*xgl*(-xb1*xb2*xt1+xb1*xb2*xt2-xb1*xt*xt1+
     . xb1*xt*xt2+xb1*xt1-xb1*xt2+xb2*xgl*xt1-xb2*xgl*xt2+xgl*xt*xt1-
     . xgl*xt*xt2-xgl*xt1+xgl*xt2)+t134(mu,mg,mt,mg)*xgl*(xb1*xgl*xt1
     . -xb1*xgl*xt2+xb1*xt*xt1-xb1*xt*xt2-xb1*xt1+xb1*xt2+xb2*xgl*xt1
     . -xb2*xgl*xt2+xb2*xt*xt1-xb2*xt*xt2-xb2*xt1+xb2*xt2-2*xgl**2*
     . xt1+2*xgl**2*xt2-2*xgl*xt*xt1+2*xgl*xt*xt2+2*xgl*xt1-2*xgl*xt2
     . )
      ans2=t134(mt1,mb1,mt,mg)*xgl*(-xb1*xb2*xt2+xb1*xb2+xb1*xgl*xt2-
     . xb1*xgl-xb2*xt*xt2+xb2*xt+xb2*xt1*xt2-xb2*xt1+xgl*xt*xt2-xgl*
     . xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mb2,mt,mg)*xgl*(-xb1*xb2*xt2+
     . xb1*xb2-xb1*xt*xt2+xb1*xt+xb1*xt1*xt2-xb1*xt1+xb2*xgl*xt2-xb2*
     . xgl+xgl*xt*xt2-xgl*xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mg,mt,mg)*
     . xgl*(xb1*xgl*xt2-xb1*xgl+xb1*xt*xt2-xb1*xt-xb1*xt1*xt2+xb1*xt1
     . +xb2*xgl*xt2-xb2*xgl+xb2*xt*xt2-xb2*xt-xb2*xt1*xt2+xb2*xt1-2*
     . xgl**2*xt2+2*xgl**2-2*xgl*xt*xt2+2*xgl*xt+2*xgl*xt1*xt2-2*xgl*
     . xt1)+t134(mt2,mb1,mt,mg)*xgl*(xb1*xb2*xt1-xb1*xb2-xb1*xgl*xt1+
     . xb1*xgl+xb2*xt*xt1-xb2*xt-xb2*xt1*xt2+xb2*xt2-xgl*xt*xt1+xgl*
     . xt+xgl*xt1*xt2-xgl*xt2)+t134(mt2,mb2,mt,mg)*xgl*(xb1*xb2*xt1-
     . xb1*xb2+xb1*xt*xt1-xb1*xt-xb1*xt1*xt2+xb1*xt2-xb2*xgl*xt1+xb2*
     . xgl-xgl*xt*xt1+xgl*xt+xgl*xt1*xt2-xgl*xt2)+ans3
      ans1=log(xb1)*log(xt1)*xb1*xt1*(-xb2*xt2+xb2+xgl*xt2-xgl)+log(
     . xb1)*log(xt2)*xb1*xt2*(xb2*xt1-xb2-xgl*xt1+xgl)+log(xb2)*log(
     . xt1)*xb2*xt1*(-xb1*xt2+xb1+xgl*xt2-xgl)+log(xb2)*log(xt2)*xb2*
     . xt2*(xb1*xt1-xb1-xgl*xt1+xgl)+log(xgl)*log(xt1)*xt1*(4*xb1*xb2
     . *xt2-4*xb1*xb2-3*xb1*xgl*xt2+3*xb1*xgl-3*xb2*xgl*xt2+3*xb2*xgl
     . +2*xgl**2*xt2-2*xgl**2)+log(xgl)*log(xt2)*xt2*(-4*xb1*xb2*xt1+
     . 4*xb1*xb2+3*xb1*xgl*xt1-3*xb1*xgl+3*xb2*xgl*xt1-3*xb2*xgl-2*
     . xgl**2*xt1+2*xgl**2)+log(xt1)**2*xt1*(-xb1*xb2*xt2+xb1*xb2+xb1
     . *xgl*xt2-xb1*xgl+xb2*xgl*xt2-xb2*xgl-xgl**2*xt2+xgl**2)+4*log(
     . xt1)*xt1*(xb1*xb2*xt2-xb1*xb2-xb1*xgl*xt2+xb1*xgl-xb2*xgl*xt2+
     . xb2*xgl+xgl**2*xt2-xgl**2)+log(xt2)**2*xt2*(xb1*xb2*xt1-xb1*
     . xb2-xb1*xgl*xt1+xb1*xgl-xb2*xgl*xt1+xb2*xgl+xgl**2*xt1-xgl**2)
     . +4*log(xt2)*xt2*(-xb1*xb2*xt1+xb1*xb2+xb1*xgl*xt1-xb1*xgl+xb2*
     . xgl*xt1-xb2*xgl-xgl**2*xt1+xgl**2)+ans2
      r52=ans1/(4*(xb1*xb2*xt1**2*xt2-xb1*xb2*xt1**2-xb1*xb2*xt1*xt2
     . **2+xb1*xb2*xt1+xb1*xb2*xt2**2-xb1*xb2*xt2-xb1*xgl*xt1**2*xt2+
     . xb1*xgl*xt1**2+xb1*xgl*xt1*xt2**2-xb1*xgl*xt1-xb1*xgl*xt2**2+
     . xb1*xgl*xt2-xb2*xgl*xt1**2*xt2+xb2*xgl*xt1**2+xb2*xgl*xt1*xt2
     . **2-xb2*xgl*xt1-xb2*xgl*xt2**2+xb2*xgl*xt2+xgl**2*xt1**2*xt2-
     . xgl**2*xt1**2-xgl**2*xt1*xt2**2+xgl**2*xt1+xgl**2*xt2**2-xgl**
     . 2*xt2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      ans5=2*((2*((14*xt2-15+14*xt1+5*xt)*xt-((xt2-3)*xt2+xt1**2+(4*
     . xt2-3)*xt1))*xgl**3+(4*xt2-3+4*xt1-8*xt)*xgl**4-((2*xt+3)*(xt-
     . xt1)**2*(xt-xt2)**2+2*xgl**5)-(2*(2*((3*xt2-2)*xt2+3*xt1**2)+(
     . 5*xt2-4)*xt1)*xt-((4*(xt2-3)*xt1-3*xt2)*xt2+(4*xt2-3)*xt1**2)-
     . 2*(24*xt2-23+24*xt1)*xt**2-10*xt**3)*xgl**2-2*((2*((3*xt2-2)*
     . xt2+3*xt1**2)+(5*xt2-4)*xt1)*xt**2+(3*((xt2-1)*xt1**2-xt2**2)+
     . (3*xt2+5)*xt1*xt2)*xt-(14*xt2-15+14*xt1)*xt**3+((xt2-3)*xt1-3*
     . xt2)*xt1*xt2+4*xt**4)*xgl)*(xt1-1)*(xt2-1)-(2*(xt+xt1)*xgl-(xt
     . -xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt-1+xgl
     . )*(xt2-2+xt1)*t134(mu,mt,mg,mg)*xgl)*(xt1-xt2)
      ans4=-2*((2*((xt-xt1)**2*(xt-xt2)**2+xgl**4-(6*xt**2+10*xt*xt1+
     . 10*xt*xt2-12*xt-xt1**2-4*xt1*xt2-xt2**2)*xgl**2-2*(xt1+xt2-xt)
     . *xgl**3+2*(((2*xt2-3)*xt2+2*xt1**2+(7*xt2-3)*xt1)*xt-(5*xt2-6+
     . 5*xt1)*xt**2-(xt1+xt2)*xt1*xt2+xt**3)*xgl)*(xt1-1)*(xt2-1)+(2*
     . (xt+xt1)*xgl-(xt-xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-
     . xgl**2)*(xt2-2+xt1)*log(xgl))*(xt1-xt2)+(2*(xt+xt1)*xgl+(xt-
     . xt1)**2*(xt1-2)-xgl**2*xt1)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2
     . )*(log(xgl)-log(xt1))*(xt2-1)**2-(2*(xt+xt1)*xgl-(xt-xt1)**2-
     . xgl**2)*(2*(xt+xt2)*xgl+(xt-xt2)**2*(xt2-2)-xgl**2*xt2)*(log(
     . xgl)-log(xt2))*(xt1-1)**2)*(log(xgl)-log(xt))*xt+ans5
      ans3=-(4*((xt-xt2**2+xt2)*(xt-xt2)**2+xgl**3-((xt2+1)*xt2+xt)*
     . xgl**2-((xt+4*xt2**2)*xt-(2*xt2-1)*xt2**2)*xgl)+(((xt+xt2)*(xt
     . -xt2)*(xt2-2)-4*xt*xt2)*xgl+(xt+2*xt2)*(xt2-2)*xgl**2-(xt-xt2)
     . **2*(xt2-2)*xt-(xt2-2)*xgl**3)*(log(xgl)-log(xt2)))*(2*(xt+xt1
     . )*xgl-(xt-xt1)**2-xgl**2)*(log(xgl)-log(xt2))*(xt1-1)**2+(4*((
     . xt-xt1**2+xt1)*(xt-xt1)**2+xgl**3-((xt1+1)*xt1+xt)*xgl**2-((xt
     . +4*xt1**2)*xt-(2*xt1-1)*xt1**2)*xgl)+(((xt+xt1)*(xt-xt1)*(xt1-
     . 2)-4*xt*xt1)*xgl+(xt+2*xt1)*(xt1-2)*xgl**2-(xt-xt1)**2*(xt1-2)
     . *xt-(xt1-2)*xgl**3)*(log(xgl)-log(xt1)))*(2*(xt+xt2)*xgl-(xt-
     . xt2)**2-xgl**2)*(log(xgl)-log(xt1))*(xt2-1)**2+ans4
      ans2=(2*(((xt+2*xt2)*xt+(xt2-4)*xt2)*xgl+(xt2+2+xt)*xgl**2-(xt+
     . xt2-2)*(xt-xt2)**2-xgl**3)*(xt1-1)**2*t134(mt2,mt,mg,mg)*xgl-(
     . 2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt2-2+xt1)*(log(xgl)+4)*(
     . xgl+xt)*(xt1-xt2)*log(xgl))*(2*(xt+xt1)*xgl-(xt-xt1)**2-xgl**2
     . )-2*(((xt+2*xt1)*xt+(xt1-4)*xt1)*xgl+(xt1+2+xt)*xgl**2-(xt+xt1
     . -2)*(xt-xt1)**2-xgl**3)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(
     . xt2-1)**2*t134(mt1,mt,mg,mg)*xgl+((2*(xt1+xt2)-xgl)*xgl**3-(xt
     . -xt1)**2*(xt-xt2)**2+(2*xt**2+6*xt*xt1+6*xt*xt2-8*xt-xt1**2-4*
     . xt1*xt2-xt2**2)*xgl**2-2*(((xt2-2)*xt2+xt1**2+2*(3*xt2-1)*xt1)
     . *xt-(3*xt2-4+3*xt1)*xt**2-(xt1+xt2)*xt1*xt2)*xgl)*(log(xgl)-
     . log(xt))**2*(xt1-xt2)*(xt1-1)*(xt2-1)*xt+ans3
      ans1=-ans2
      r82=ans1/(4*((xt-xt1)**2+xgl**2-2*(xt+xt1)*xgl)*((xt-xt2)**2+
     . xgl**2-2*(xt+xt2)*xgl)*(xt1-xt2)*(xt1-1)**2*(xt2-1)**2)

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      ans6=-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(((xt1-xt2)*xt1+
     . 4*xt**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+4*log(xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2-4*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-
     . xt1)**2+xgl**2)))*xt1+(xt1-xt2)*log(-(xt-xt1+xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1)*(xt-xt1+xgl)*xt2
      ans5=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(((xt1-xt2)*xt2-
     . 4*xt**2)*(log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . )-4*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*
     . xt**2+4*log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2)))*xt2+(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt2)*(xt-xt2+xgl)*xt1+ans6
      ans4=4*(xgl+xt-xt1)*(xgl-xt-xt1)*log(xt-xt1-xgl+sqrt(-2*(xt+xt1
     . )*xgl+(xt-xt1)**2+xgl**2))*xt**2*xt2+4*(xgl+xt-xt1)*(xgl-xt-
     . xt1)*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl
     . +sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-
     . xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*
     . xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(
     . xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(xt1-xt2)*log(-(xt-
     . xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1*xt2+
     . ans5
      ans3=(2*((xt1+xt2-2*xt-2*xgl)*(xt1-xt2)+(xt1-xt2-8*xt)*log(xt2)
     . *xt2+(xt1-xt2+8*xt)*log(xt1)*xt1-6*(xt1-xt2)*log(xt)*xt-(xt1+
     . xt2+2*xt)*(xt1-xt2)*log(xgl))*xt*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)
     . *(4*xt**2-xt1*xt2+xt2**2)*log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+
     . (xt-xt2)**2+xgl**2))+(xgl+xt-xt2)*(xgl-xt-xt2)*(4*xt**2-xt1*
     . xt2+xt2**2)*log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2)))*xt1-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*
     . xt2)*log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt2-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*xt2)*log(xt+
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*xt2-4*(xgl+
     . xt-xt2)*(xgl-xt-xt2)*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2))*xt**2*xt1-4*(xgl+xt-xt2)*(xgl-xt-xt2)*log(xt-
     . xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*xt**2*xt1+
     . ans4
      ans7=((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)
      ans2=ans3*ans7
      ans1=-ans2
      rat2=ans1/(16*(xt1-xt2)**2*(xt1-1)*(xt2-1)*xt**2*xt1*xt2)

      ans5=((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . )))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-xt2)*log(xgl)
      ans4=sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(-(xt-xt2+xgl+sqrt(-2*
     . (xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))+log(-(xt-xt2+xgl-sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))))*(xt-xt2+xgl)*((xt1-1)*log(
     . xt2)*xt2-(xt2-1)*log(xt1)*xt1)+2*((log(xgl)-log(xt2))**2*(2*xt
     . -xt2)*(xt1-1)*xt2-2*(xt1-xt2)*log(xgl)**2*xt-(log(xgl)-log(xt1
     . ))**2*(2*xt-xt1)*(xt2-1)*xt1)*xt+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)
     . **2+xgl**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+
     . xgl**2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**
     . 2))-log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+log(-(xt-xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))))
     . *(xt-xt1+xgl)*((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)+ans5
      ans3=-ans4
      ans2=-(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)
     . )+log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt2-1)-2*((log(xgl)-log(xt2))*xt2-log(xgl))*(xt1-xt2)*xt)*(log
     . (xgl)-log(xt1))*xt1+ans3
      ans1=(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-1)-2*(xt1-xt2)*log(xgl)*xt)*(log(xgl)-log(xt2))*xt2+ans2
      rlt2=ans1/(8*(xt1-xt2)*(xt1-1)*(xt2-1)*xt**2)

      ans3=(2*(((log(xgl)-log(xt2))*(xt1-1)**2*xt2-(xt1-xt2)**2*log(
     . xgl)+(xt2+1-7*xt)*(xt2-1)*xt1)*xt1-(7*xt1**2-3*xt2-2*(xt2+1)*
     . xt1)*(xt2-1)*xgl)*xt1+(2*((2*(xt2+1)*xt1+3*xt2)*xt-2*xt1**3)*
     . xt1-(log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)*(xt1**2-xt2))*(xt2-1))*(log(xgl)-log(xt1))
      ans2=sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(log(xt+xt1-xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+sqrt(-2*(
     . xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt-xt1-xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*((xt1**2-xt2)*(xt2-1)*log(xt1)-
     . (xt1-1)**2*log(xt2)*xt2-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt1+xgl
     . )+((log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)+2*(5*xt+xt1+5*xgl)*xt1)*((log(xgl)-log(xt2))*(
     . xt1-1)**2*xt2-(xt1-xt2)**2*log(xgl))+ans3
      ans1=2*(((log(xgl)-log(xt2))**2*(xt1-1)**2*xt2-(xt1-xt2)**2*log
     . (xgl)**2)*(xgl+xt)+(xt2+1-5*xt)*(xt2-1)*xt1**2)*xt1-(log(xt-
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-
     . xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl
     . )*(xt1-xt2)*(xt1-1)*(xt2-1)-2*(xt+xt1+xgl)*(log(xgl)-log(xt1))
     . **2*(xt1**2-xt2)*(xt2-1)*xt1-2*(5*xt*xt2+xt1**3+5*(xt1-xt2)*(
     . xt1-1)*xgl-(5*(xt2+1)*xt-xt2)*xt1)*(xt2-1)*xt1-2*((log(xgl)-
     . log(xt1))*(xt1**2-xt2)*(xt2-1)-(log(xgl)-log(xt2))*(xt1-1)**2*
     . xt2+((xt1-xt2)*log(xgl)+(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt1+ans2
      rmt12=ans1/(4*(xt1-xt2)**2*(xt1-1)**2*(xt2-1)*xt1)

      ans3=((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . -log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(log(xgl)-log(xt2))*xt2**2+2*(5*xt+xt2+5*xgl)
     . *xt2)*(log(xgl)-log(xt1))*(xt2-1)**2*xt1+(2*(((2*xt2+3)*xt+xt2
     . **2)*xt1**2-(xt1-xt2)**2*log(xgl)*xt2-((7*xt2**2+3)*xt+2*xt2**
     . 3)*xt1)*xt2+(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2))-log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2))+log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)*(xt1-xt2**2)*(xt1-1)+2*(((7*xt-1)*xt2-2*
     . (xt-xt2**2))*xt2-((7*xt2-2)*xt2-(2*xt2+3)*xt1)*(xt1-1)*xgl)*
     . xt2)*(log(xgl)-log(xt2))
      ans2=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(xt-xt2-xgl+sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt-xt2-xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*((xt1-xt2**2)*(xt1-1)*log(xt2)+
     . (xt2-1)**2*log(xt1)*xt1-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt2+xgl
     . )-((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(5*xt+xt2+5*xgl)*xt2)*(xt1-xt2)**2*log(xgl)+
     . ans3
      ans1=(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)+2*((log(xgl)-log(xt1)
     . )**2*(xt2-1)**2*xt1-(xt1-xt2)**2*log(xgl)**2)*(xgl+xt)*xt2+2*(
     . xt+xt2+xgl)*(log(xgl)-log(xt2))**2*(xt1-xt2**2)*(xt1-1)*xt2+2*
     . (5*xt+xt2+5*xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)*xt2+2*((log(xgl)-
     . log(xt1))*(xt2-1)**2*xt1+(log(xgl)-log(xt2))*(xt1-xt2**2)*(xt1
     . -1)-((xt1-xt2)*log(xgl)-(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt2+ans2
      rmt22=ans1/(4*(xt1-xt2)**2*(xt1-1)*(xt2-1)**2*xt2)

      relw = r22+r42+r52+r72+r82+r92+rlt2+rmt12+rmt22
      bo  = fi(amst1**2,amst2**2,amu**2)
      relw = cf*relw/bo

      anomalous =-cf
      finscale = 2*dlog(scale**2/amg**2)
      felw_hdec = dreal(relw)*fnorm + anomalous + finscale

      return
	  
      end

! ======================================================================================================
! ====================== Lay down old test for the main subroutine code here ===========================
! ======================================================================================================

c	integer a1, a2, a3, a4			! To test the stability of t134_hdec
c
c	do a1=1,5					
c		do a2=1,5
c			do a3=1,5
c				do a4=1,5
c				t134testC = t134_hdec(a1*100.0d0, a2*100.0d0, a3*100.0d0, a4*100.0d0)
c				write(*,*)'t134C Test lieferte: ', t134testC, ' für a1, a2, a3, a4: ', a1,a2,a3,a4
c				enddo
c			enddo
c		enddo
c	enddo

c	common/poly_hdec/nber			! Needed for the Dilog of HDECAY instead of cspen(x) from DM@NLO
c	nber = 18
c	call bernini_hdec(nber)																					
			  
! ======================================================================================================
! ========= Code below is needed for the Dilog of HDECAY instead of cspen(x) from DM@NLO ===============								
! ======================================================================================================


c      double precision function sp_hdec(x)									! From hdecay.f l.10576
	  ! Real part of the Dilogarithm
c	  
c      implicit double precision (a-h,o-z)
c      complex*16 cx,li2_hdec
c      cx = dcmplx(x,0.d0)
c      sp_hdec = dreal(li2_hdec(cx))
c      return
c	  
c      end
 
! ------------------------------------------------------------------------------------------------------
 
 
c      complex*16 function li2_hdec(x)										! From hdecay.f l.10585
	  ! Complex Dilogarithm
	  
c      implicit double precision (a-h,o-z)
c      complex*16 x,y,cli2_hdec
c      common/const_hdec/zeta2,zeta3
c      zero=1.d-16
c      xr=dreal(x)
c      xi=dimag(x)
c      r2=xr*xr+xi*xi
c      li2_hdec=0
c      if(r2.le.zero)then
c        li2_hdec=x
c        return
c      endif
c      rr=xr/r2
c      if(r2.eq.1.d0.and.xi.eq.0.d0)then
c        if(xr.eq.1.d0)then
c          li2_hdec=dcmplx(zeta2)
c        else
c          li2_hdec=-dcmplx(zeta2/2.d0)
c        endif
c        return
c      elseif(r2.gt.1.d0.and.rr.gt.0.5d0)then
c        y=(x-1.d0)/x
c        li2_hdec=cli2_hdec(y)+zeta2-cdlog(x)*cdlog(1.d0-x)
c     .          +0.5d0*cdlog(x)**2
c        return
c      elseif(r2.gt.1.d0.and.rr.le.0.5d0)then
c        y=1.d0/x
c        li2_hdec=-cli2_hdec(y)-zeta2-0.5d0*cdlog(-x)**2
c        return
c      elseif(r2.le.1.d0.and.xr.gt.0.5d0)then
c        y=1.d0-x
c        li2_hdec=-cli2_hdec(y)+zeta2-cdlog(x)*cdlog(1.d0-x)
c       return
c      elseif(r2.le.1.d0.and.xr.le.0.5d0)then
c        y=x
c        li2_hdec=cli2_hdec(y)
c        return
c      endif
c
c      end
	  
! ------------------------------------------------------------------------------------------------------
 
 
c      complex*16 function cli2_hdec(x)										! From hdecay.f l.10627
	  ! Taylor Expansion of the Dilogarithm
	  
c      implicit double precision (a-h,o-z)
c      complex*16 x,z
c      common/bernoulli_hdec/b2(18),b12(18),b3(18)
c      common/poly_hdec/nber
c      n=nber-1
c      z=-cdlog(1.d0-x)
c      cli2_hdec=b2(nber)
c      do 111 i=n,1,-1
c        cli2_hdec=z*cli2_hdec+b2(i)
c111   continue
c      cli2_hdec=z**2*cli2_hdec+z
c      return
c	  
c      end
 
! ------------------------------------------------------------------------------------------------------
 
 
c      double precision function factrl_hdec(n)								! From hdecay.f l.10643
	  ! Factorial function
	  
c      implicit double precision (a-h,o-z)
c      factrl_hdec=1.d0
c      if(n.eq.0)return
c      do 999 i=1,n
c        factrl_hdec=factrl_hdec*dfloat(i)
c999   continue
c      return
c	  
c      end
 
! ------------------------------------------------------------------------------------------------------
 
 
c      subroutine bernini_hdec(n)											! From hdecay.f l.10654
	  ! Initialization of Bernoulli numbers
	  
c      implicit double precision (a-h,o-z)
c      dimension b(18),pb(19)
c      common/bernoulli_hdec/b2(18),b12(18),b3(18)
c      common/const_hdec/zeta2,zeta3
c      common/poly_hdec/nber
c 
c      nber=n
c      pi=4.d0*datan(1.d0)
c 
c      b(1)=-1.d0/2.d0
c      b(2)=1.d0/6.d0
c      b(3)=0.d0
c      b(4)=-1.d0/30.d0
c      b(5)=0.d0
c      b(6)=1.d0/42.d0
c      b(7)=0.d0
c      b(8)=-1.d0/30.d0
c      b(9)=0.d0
c      b(10)=5.d0/66.d0
c      b(11)=0.d0
c      b(12)=-691.d0/2730.d0
c      b(13)=0.d0
c      b(14)=7.d0/6.d0
c      b(15)=0.d0
c      b(16)=-3617.d0/510.d0
c      b(17)=0.d0
c      b(18)=43867.d0/798.d0
c      zeta2=pi**2/6.d0
c      zeta3=1.202056903159594d0
c 
c      do 995 i=1,18
c        b2(i)=b(i)/factrl_hdec(i+1)
c        b12(i)=dfloat(i+1)/factrl_hdec(i+2)*b(i)/2.d0
c        pb(i+1)=b(i)
c        b3(i)=0.d0
c995   continue
c      pb(1)=1.d0
c      do 996 i=1,18
c      do 996 j=0,i
c        b3(i)=b3(i)+pb(j+1)*pb(i-j+1)/factrl_hdec(i-j)/factrl_hdec(j+1)
c     .                                            /dfloat(i+1)
c996   continue
c 
c      return
c	  
c      end
	
	