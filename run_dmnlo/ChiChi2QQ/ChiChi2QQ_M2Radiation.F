************************************************

	double precision function ChiChi2QQ_RadiationS()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j,imin,imax,jmin,jmax

	double precision temp1,temp2
	complex*16 DenSmZ, DenSmZc, DenSmH(6)
	double precision DenSmZ1, DenSmH1(6)
	
***** widths of Z-boson and Higgs included in propagators	

	if (icha.eq.1) then

	DenSmZ = (sqrtS**2 - MW2) + cI*MW*wW
        DenSmZc = (sqrtS**2 - MW2) - cI*MW*wW
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MW2)
		
	else
		
	DenSmZ = (sqrtS**2 - MZ2) + cI*MZ*wZ
        DenSmZc = (sqrtS**2 - MZ2) - cI*MZ*wZ
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MZ2)

	endif


	if (icha.eq.1) then
c only charged Higgs exchange		
	  imin = 5
	  imax = 6 
	  jmin = 5
	  jmax = 6
	else
c exchange of neutral Higgses 		
	  imin = 1
	  imax = 4 
	  jmin = 1
	  jmax = 4	  	
	endif

	do i = imin, imax
	  DenSmH(i) = (sqrtS**2 - MHk(i)**2) + cI*MHk(i)*wHk(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk(i)**2)
	enddo


! Photon contributions only in cha-cha case
! cha cha -> b B : no t-channel
! cha cha -> t T : no u-channel
! Use couplings E/F, G/H to switch on/off t-/u-channels	


        ChiChi2QQ_RadiationS = 0d0 
        
** Photon exchange in Chargino annihilation (m1=m2, m3=m4)

	if (icha.eq.2) then
		
** Photon 

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS +  1/(4.*pi)*  (-128*ApL**2*BpL**2*CF*gs**2*NC*Pi*(cosxi*(-1 + 4*muu1**2)*muu3**2*(-1 + x1)
     &       *Sqrt(-4*muu3**2 + x1**2)*(-1 + x2)*Sqrt(-4*muu3**2 + x2**2) - 
     &      (1 + 2*muu1**2)*(-1 + x1)*(-1 + x2)*(x1**2 + x2**2) + 4*muu3**4*(5 + x1**2 - 5*x2 + x2**2 + x1*(-5 + 3*x2) + 2*muu1**2*(2 - 2*x1 + x1**2 - 2*x2 + x2**2)) + 
     &      muu3**2*(14 - 20*x2 + 8*x2**2 + x1*(-20 + 23*x2 - 7*x2**2) + x1**2*(8 - 7*x2 + x2**2) + 4*muu1**2*((-2 + x2)**2 + x1**2*(1 + x2 - x2**2) + x1*(-4 + x2 + x2**2)))))/
     &  (3.*sqrtS**2*(-1 + x1)**2*(-1 + x2)**2)

        if (cutoff.lt.zerocheck) then

** Photon dipole contribution 312                              

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - 1/(4.*pi) * (128*ApL**2*BpL**2*CF*GS**2*(1 + 2*muu1**2)*(1 + 2*muu3**2)*NC*Pi*(4*muu3**2 - x2**2)*
     &    ((-2 + 4*muu3**2)/(-2 + x1 + x2) + (Sqrt((-1 + 4*muu3**2)/(4*muu3**2 - x2**2))*(4*muu3**4 + muu3**2*(4 - 6*x2) + (-1 + x2)*(-1 + x1 + 2*x2)))/((-1 + 2*muu3**2)*(-1 + x2))))/
     &  (3.*sqrtS**2*(-1 + x2)*(-4*muu3**2 + x2**2))
    

** Photon dipole contribution 321

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - 1/(4.*pi)* (128*ApL**2*BpL**2*CF*GS**2*(1 + 2*muu1**2)*(1 + 2*muu3**2)*NC*Pi*(4*muu3**2 - x1**2)*
     &    ((-2 + 4*muu3**2)/(-2 + x1 + x2) + (Sqrt((-1 + 4*muu3**2)/(4*muu3**2 - x1**2))*(4*muu3**4 + muu3**2*(4 - 6*x1) + (-1 + x1)*(-1 + 2*x1 + x2)))/((-1 + 2*muu3**2)*(-1 + x1))))/
     &  (3.*sqrtS**2*(-1 + x1)*(-4*muu3**2 + x1**2))

        endif
    
** G-Z interference

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS +  1/(4.*pi)*2d0*dreal( (-4*ApL*(AL + AR)*BpL*(BL + BR)*CF*gs**2*NC*Pi*(-((-1 + 2*cosxi**2)
     &      *(-1 + 4*muu1**2)*(4*muu3**2 - x1**2)*(1 - x1 + 2*muu3**2*(-1 + x2))*(-1 + x2)) + 
     &      (-1 + x1)*(-1 + x2)*((-7 - 20*muu1**2 + cosxi**2*(-2 + 8*muu1**2))*x1**2 - 8*(1 + 2*muu1**2)*x2**2) + 
     &      2*muu3**2*(54 + 72*muu1**2 - 78*x1 - 72*muu1**2*x1 + 31*x1**2 + 20*muu1**2*x1**2 - 2*cosxi**2*(-1 + 4*muu1**2)*(-2 + 2*x1 + x1**2*(-1 + x2))*(-1 + x2) - 78*x2 - 72*muu1**2*x2 + 
     &         90*x1*x2 + 24*muu1**2*x1*x2 - 26*x1**2*x2 + 8*muu1**2*x1**2*x2 + 32*x2**2 + 16*muu1**2*x2**2 - 28*x1*x2**2 + 16*muu1**2*x1*x2**2 + 3*x1**2*x2**2 - 12*muu1**2*x1**2*x2**2 + 
     &         4*cosxi*(-1 + 4*muu1**2)*(-1 + x1)*Sqrt(-4*muu3**2 + x1**2)*(-1 + x2)*Sqrt(-4*muu3**2 + x2**2)) + 
     &      8*muu3**4*(21 - 20*x1 + 4*x1**2 + 2*cosxi**2*(-1 + 4*muu1**2)*(-1 + x2)**2 - 22*x2 + 12*x1*x2 + 5*x2**2 + 4*muu1**2*(3 - 4*x1 + 2*x1**2 - 2*x2 + x2**2))))/
     &  dreal(3.*DenSmZ*(-1 + x1)**2*(-1 + x2)**2) )  


        if (cutoff.lt.zerocheck) then

** G-Z-Boson dipole contribution 312                             
                                                                
	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - 1/(4.*pi)*2d0*dreal(  (-32*ApL*(AL + AR)*BpL*(BL + BR)*CF*GS**2*(1 + 2*muu1**2)*(1 + 2*muu3**2)*NC*Pi*
     &    ((-1 + x2)*(2 + 2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + x1**2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + 3*x1*(-1 + x2)*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) - 
     &         5*x2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + 2*x2**2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2))) + 
     &      4*muu3**4*(-2 - 2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + x1*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + x2*(2 + Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)))) - 
     &      2*muu3**2*(-4 + 4*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) - 2*x1*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + 3*x2**2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + 
     &         x2*(4 - 8*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2)) + 3*x1*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x2**2))))))/dreal(3.*DenSmZ*(-1 + 2*muu3**2)*(-1 + x2)**2*(-2 + x1 + x2))  )   

** G-Z-Boson dipole contribution 321                             

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - 1/(4.*pi)* 2d0*dreal( (-32*ApL*(AL + AR)*BpL*(BL + BR)*CF*GS**2*(1 + 2*muu1**2)*(1 + 2*muu3**2)*NC*Pi*
     &    (4*muu3**4*(-2 - 2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2)) + x1*(2 + Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2))) + Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2))*x2) + 
     &      (-1 + x1)*(2 + 2*Sqrt((-1 + 4*muu3**2)/(4*muu3**2 - x1**2)) + 2*x1**2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2)) - 3*Sqrt((-1 + 4*muu3**2)/(4*muu3**2 - x1**2))*x2 + 
     &         Sqrt((-1 + 4*muu3**2)/(4*muu3**2 - x1**2))*x2**2 + x1*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2))*(-5 + 3*x2)) - 
     &      2*muu3**2*(-4 + 4*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2)) + 3*x1**2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2)) - 2*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2))*x2 + 
     &         x1*(4 - 8*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2)) + 3*Sqrt((1 - 4*muu3**2)/(-4*muu3**2 + x1**2))*x2))))/dreal(3.*DenSmZ*(-1 + 2*muu3**2)*(-1 + x1)**2*(-2 + x1 + x2)) )    

        endif



	endif

** Z/W-boson bremsstrahlung 11+22+21+12

	ChiChi2QQ_RadiationS= ChiChi2QQ_RadiationS +  (8*CF*gs**2*NC*(((2*AR*muu1*muu2*sqrtS**4*(4*BL*BR*muu3*muu4*(2*muu3**2 + 2*muu4**2 - x1 - x2) + 
     &               BL**2*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2)) + 
     &               BR**2*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))) + 
     &            AL*(-2*BL*BR*muu3*muu4*(8*p1k1**2 + 8*p1k2**2 - 4*p1k2*sqrtS**2*(2*muu1**2 - 2*muu2**2 + x1 + x2) + 
     &                  4*p1k1*(4*p1k2 - sqrtS**2*(2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(2*muu1**4 + (-1 + muu2**2)*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 3*x2) + 
     &                     muu1**2*(-4*muu2**2 + 2*muu3**2 + 2*muu4**2 + x1 + x2))) + 
     &               BL**2*(8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 
     &                  2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 2*muu2**2*muu3**2*p1k2*sqrtS**2 - 
     &                  4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + 
     &                  sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + 
     &                  muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + 
     &                  muu4**2*(8*p1k1**2 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1) + 
     &                     sqrtS**2*((2 - 2*muu1**2 + 2*muu2**2)*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)) - 
     &                  4*p1k1**2*x2 - 2*p1k2*sqrtS**2*x2 + 2*muu1**2*p1k2*sqrtS**2*x2 - 2*muu2**2*p1k2*sqrtS**2*x2 - 
     &                  4*muu3**2*p1k2*sqrtS**2*x2 + 4*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 - muu1**2*sqrtS**4*x1*x2 + 
     &                  muu2**2*sqrtS**4*x1*x2 - 2*p1k1*
     &                   (2*p1k2*(-2 + x1 + x2) + sqrtS**2*
     &                      (1 + muu3**2 - x1 + muu1**2*(1 + muu3**2 - x1 - x2) - x2 - x1*x2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2)))) + 
     &               BR**2*(-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*
     &                   (muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu1**2*(1 + muu3**2 - x1 - x2) + 3*x2 - 2*x1*x2 - x2**2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2))))))*dconjg(AL) + 
     &         (2*AL*muu1*muu2*sqrtS**4*(4*BL*BR*muu3*muu4*(2*muu3**2 + 2*muu4**2 - x1 - x2) + 
     &               BL**2*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2)) + 
     &               BR**2*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))) + 
     &            AR*(-2*BL*BR*muu3*muu4*(8*p1k1**2 + 8*p1k2**2 - 4*p1k2*sqrtS**2*(2*muu1**2 - 2*muu2**2 + x1 + x2) + 
     &                  4*p1k1*(4*p1k2 - sqrtS**2*(2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(2*muu1**4 + (-1 + muu2**2)*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 3*x2) + 
     &                     muu1**2*(-4*muu2**2 + 2*muu3**2 + 2*muu4**2 + x1 + x2))) + 
     &               BR**2*(8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 
     &                  2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 2*muu2**2*muu3**2*p1k2*sqrtS**2 - 
     &                  4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + 
     &                  sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + 
     &                  muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + 
     &                  muu4**2*(8*p1k1**2 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1) + 
     &                     sqrtS**2*((2 - 2*muu1**2 + 2*muu2**2)*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)) - 
     &                  4*p1k1**2*x2 - 2*p1k2*sqrtS**2*x2 + 2*muu1**2*p1k2*sqrtS**2*x2 - 2*muu2**2*p1k2*sqrtS**2*x2 - 
     &                  4*muu3**2*p1k2*sqrtS**2*x2 + 4*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 - muu1**2*sqrtS**4*x1*x2 + 
     &                  muu2**2*sqrtS**4*x1*x2 - 2*p1k1*
     &                   (2*p1k2*(-2 + x1 + x2) + sqrtS**2*
     &                      (1 + muu3**2 - x1 + muu1**2*(1 + muu3**2 - x1 - x2) - x2 - x1*x2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2)))) + 
     &               BL**2*(-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*
     &                   (muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu1**2*(1 + muu3**2 - x1 - x2) + 3*x2 - 2*x1*x2 - x2**2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2))))))*dconjg(AR))/
     &       ((1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)) - 
     &      (AL*((BR*(-(BL*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 + muu4**2 - x2)) + 
     &                  BR*(2*p1k2 - sqrtS**2*x2)*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                     (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))) + 
     &               BL*(-(BR*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 + muu4**2 - x2)) + 
     &                  2*BL*p1k2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 2*p1k2*(-1 + 3*muu3**2 - muu4**2 + x2) - 
     &                     sqrtS**2*(1 - 3*muu3**2 + muu4**2 - x1 + muu3**2*x1 - muu4**2*x1 + 
     &                        muu2**2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*x2 + 3*muu3**2*x2 - muu4**2*x2 + x1*x2 + x2**2 + 
     &                        muu1**2*(-1 + 3*muu3**2 - muu4**2 + x2)))))*dconjg(AL) + 
     &            muu1*muu2*sqrtS**4*(BL**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + 
     &                  (-1 + muu4**2 + x1)*(1 + muu4**2 - x2)) + 
     &               BR**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2)) + 
     &               8*BL*BR*muu3*muu4*(1 + muu3**2 + muu4**2 - x2))*dconjg(AR)) + 
     &         AR*(muu1*muu2*sqrtS**4*(BL**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + 
     &                  (-1 + muu4**2 + x1)*(1 + muu4**2 - x2)) + 
     &               BR**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2)) + 
     &               8*BL*BR*muu3*muu4*(1 + muu3**2 + muu4**2 - x2))*dconjg(AL) + 
     &            (BL*(-(BR*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 + muu4**2 - x2)) + 
     &                  BL*(2*p1k2 - sqrtS**2*x2)*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                     (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))) + 
     &               BR*(-(BL*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 + muu4**2 - x2)) + 
     &                  2*BR*p1k2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 2*p1k2*(-1 + 3*muu3**2 - muu4**2 + x2) - 
     &                     sqrtS**2*(1 - 3*muu3**2 + muu4**2 - x1 + muu3**2*x1 - muu4**2*x1 + 
     &                        muu2**2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*x2 + 3*muu3**2*x2 - muu4**2*x2 + x1*x2 + x2**2 + 
     &                        muu1**2*(-1 + 3*muu3**2 - muu4**2 + x2)))))*dconjg(AR)))/(-1 + muu3**2 - muu4**2 + x2)**2 - 
     &      (-4*(BL**2 + BR**2)*p1k1**2*(1 + muu3**2 - 3*muu4**2 - x1)*(AL*dconjg(AL) + AR*dconjg(AR)) + 
     &         2*p1k1*(-2*(BL**2 + BR**2)*p1k2*(1 + muu3**2 - muu4**2 - x1)*(AL*dconjg(AL) + AR*dconjg(AR)) + 
     &            sqrtS**2*(AL*(BL**2*(1 + muu3**2 - 3*muu4**2 - x1)*(1 + muu1**2 - muu2**2 + x1) + 
     &                  BR**2*(-1 - muu3**2 + 3*muu4**2 + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1) + 2*x1 + muu3**2*x1 - 
     &                     3*muu4**2*x1 - x1**2 + muu2**2*(-1 - muu3**2 + 3*muu4**2 + x1) + x2 + muu3**2*x2 - muu4**2*x2 - 
     &                     x1*x2))*dconjg(AL) + AR*
     &                (BR**2*(1 + muu3**2 - 3*muu4**2 - x1)*(1 + muu1**2 - muu2**2 + x1) + 
     &                  BL**2*(-1 - muu3**2 + 3*muu4**2 + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1) + 2*x1 + muu3**2*x1 - 
     &                     3*muu4**2*x1 - x1**2 + muu2**2*(-1 - muu3**2 + 3*muu4**2 + x1) + x2 + muu3**2*x2 - muu4**2*x2 - 
     &                     x1*x2))*dconjg(AR))) + sqrtS**2*
     &          (2*p1k2*(1 + muu3**2 - muu4**2 - x1)*x1*(AL*BL**2*dconjg(AL) + AR*BR**2*dconjg(AR)) + 
     &            sqrtS**2*((-(AL*BL*(2*BR*(-1 + muu1**2 + muu2**2)*muu3*muu4*(1 + muu3**2 + muu4**2 - x1) + 
     &                       BL*(1 + muu1**2 - muu2**2)*(1 + muu3**2 - 3*muu4**2 - x1)*x1)) + 
     &                  AR*muu1*muu2*(8*BL*BR*muu3*muu4*(1 + muu3**2 + muu4**2 - x1) + 
     &                     BL**2*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*(3*x1 - x2) + x2 - x1*x2 + 
     &                        muu3**2*(-6*muu4**2 - x1 + x2)) + 
     &                     BR**2*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*(3*x1 - x2) + x2 - x1*x2 + 
     &                        muu3**2*(-6*muu4**2 - x1 + x2))))*dconjg(AL) + 
     &               (-(AR*BR*(2*BL*(-1 + muu1**2 + muu2**2)*muu3*muu4*(1 + muu3**2 + muu4**2 - x1) + 
     &                       BR*(1 + muu1**2 - muu2**2)*(1 + muu3**2 - 3*muu4**2 - x1)*x1)) + 
     &                  AL*muu1*muu2*(8*BL*BR*muu3*muu4*(1 + muu3**2 + muu4**2 - x1) + 
     &                     BL**2*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*(3*x1 - x2) + x2 - x1*x2 + 
     &                        muu3**2*(-6*muu4**2 - x1 + x2)) + 
     &                     BR**2*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*(3*x1 - x2) + x2 - x1*x2 + 
     &                        muu3**2*(-6*muu4**2 - x1 + x2))))*dconjg(AR))))/(-1 - muu3**2 + muu4**2 + x1)**2))/
     &  (DenSmZ*sqrtS**2*dconjg(DenSmZ))

        if (cutoff.lt.zerocheck) then
       
** Zboson dipole contribution 312                                                                               

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS -  (8*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (-(((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &             ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &             (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &             ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &              (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &           (-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4) + 
     &             sqrtS**2*((1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2 + 
     &                (-1 + muu1**2 - muu2**2)*(-1 + muu3**2 - muu4**2)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &           (AL*BR**2*dconjg(AL) + AR*BL**2*dconjg(AR)))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) - 
     &      ((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. + 
     &           ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &           (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &           ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &            (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &         (4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4) + 
     &           sqrtS**2*(-((1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2) - 
     &              (-1 + muu1**2 - muu2**2)*(1 + muu3**2 - muu4**2)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &         (AL*BL**2*dconjg(AL) + AR*BR**2*dconjg(AR)))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &      sqrtS**4*((BL**2 + BR**2)*muu1*muu2*(-1 + muu3**2)*(AR*dconjg(AL) + AL*dconjg(AR)) + 
     &         (BL**2 + BR**2)*muu1*muu2*muu4**2*(AR*dconjg(AL) + AL*dconjg(AR)) + 
     &         2*BL*BR*muu3*muu4*((-4*AR*muu1*muu2 + AL*(-1 + muu1**2 + muu2**2))*dconjg(AL) + 
     &            (-4*AL*muu1*muu2 + AR*(-1 + muu1**2 + muu2**2))*dconjg(AR)))))/
     &  (DenSmZ*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmZ))

** Zboson dipole contribution 321

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - (-8*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (3*muu3**4 - muu4**4 + muu3**2*(4 + 2*muu4**2 - 5*x1 - x2) + muu4**2*(-x1 + x2) + (-1 + x1)*(-1 + 2*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (-((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &            (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &               x1*((-1 + muu4**2)*sqrtS**4*x1 + Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                   Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &          (-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &            sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - 
     &               (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &               x1*((-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**4*x1 + 
     &                  (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                   Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*(AL*BR**2*dconjg(AL) + AR*BL**2*dconjg(AR)))/
     &       (4.*sqrtS**8*(-4*muu3**2 + x1**2)**2) - ((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &            Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &              x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                  Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &         (-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - 
     &              (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &              x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &                 (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                  Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*(AL*BL**2*dconjg(AL) + AR*BR**2*dconjg(AR)))/
     &       (4.*sqrtS**8*(-4*muu3**2 + x1**2)**2) + sqrtS**4*
     &       ((BL**2 + BR**2)*muu1*muu2*(-1 + muu3**2)*(AR*dconjg(AL) + AL*dconjg(AR)) + 
     &         (BL**2 + BR**2)*muu1*muu2*muu4**2*(AR*dconjg(AL) + AL*dconjg(AR)) + 
     &         2*BL*BR*muu3*muu4*((-4*AR*muu1*muu2 + AL*(-1 + muu1**2 + muu2**2))*dconjg(AL) + 
     &            (-4*AL*muu1*muu2 + AR*(-1 + muu1**2 + muu2**2))*dconjg(AR)))))/
     &  (DenSmZ*sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmZ))

        endif

	do i = imin, imax

** Z/W-Higgs interference                                                      

	  ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS + 1/(4.*pi)* 2*dreal(  (-8*CF*gs**2*NC*Pi*sqrtS**2*(CR(i)*((muu1 - muu1**3 + muu1*muu2**2)
     &      *dconjg(AL) + muu2*(-1 - muu1**2 + muu2**2)*dconjg(AR)) + 
     &      CL(i)*(muu2*(-1 - muu1**2 + muu2**2)*dconjg(AL) + muu1*(1 - muu1**2 + muu2**2)*dconjg(AR)))*
     &    ((BR*muu4*(2*muu3**6 + muu3**4*(10 - 6*muu4**2 + x1**2 + 4*x1*(-2 + x2) - 8*x2 + 3*x2**2) + 
     &            muu3**2*(6 + 6*muu4**4 - x1**3 - 3*x1**2*(-2 + x2) - 4*x2 + x2**3 - x1*(12 - 8*x2 + x2**2) - 4*muu4**2*(3 + x1*(-2 + x2) - 2*x2 + x2**2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(-4 + x1**2 + 2*x2 - x2**2) - (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) - 2*x2 + x2**2))) + 
     &         BL*muu3*(2*muu3**6 - muu3**4*(2 + 6*muu4**2 + x1**2 - x2**2) + 
     &            muu3**2*(-2 + 6*muu4**4 + x1**3 + 4*muu4**2*(3 + x1**2 + x1*(-2 + x2) - 2*x2) - 4*x2 + 4*x2**2 - x2**3 - x1**2*(2 + x2) + x1*(4 + 4*x2 - 3*x2**2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(8 + 3*x1**2 + 4*x1*(-2 + x2) - 6*x2 + x2**2) + (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) - 2*x2 + x2**2))))*DL(i) + 
     &      (BL*muu4*(2*muu3**6 + muu3**4*(10 - 6*muu4**2 + x1**2 + 4*x1*(-2 + x2) - 8*x2 + 3*x2**2) + 
     &            muu3**2*(6 + 6*muu4**4 - x1**3 - 3*x1**2*(-2 + x2) - 4*x2 + x2**3 - x1*(12 - 8*x2 + x2**2) - 4*muu4**2*(3 + x1*(-2 + x2) - 2*x2 + x2**2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(-4 + x1**2 + 2*x2 - x2**2) - (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) - 2*x2 + x2**2))) + 
     &         BR*muu3*(2*muu3**6 - muu3**4*(2 + 6*muu4**2 + x1**2 - x2**2) + 
     &            muu3**2*(-2 + 6*muu4**4 + x1**3 + 4*muu4**2*(3 + x1**2 + x1*(-2 + x2) - 2*x2) - 4*x2 + 4*x2**2 - x2**3 - x1**2*(2 + x2) + x1*(4 + 4*x2 - 3*x2**2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(8 + 3*x1**2 + 4*x1*(-2 + x2) - 6*x2 + x2**2) + (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) - 2*x2 + x2**2))))*DR(i)))/
     &  dreal((-1 - muu3**2 + muu4**2 + x1)**2*(-1 + muu3**2 - muu4**2 + x2)**2*dconjg(DenSmZ)*DenSmH(i))   )    

        if (cutoff.lt.zerocheck) then
   	 
** Z-Higgs dipole contribution 312                                             

       ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS -1/(4.*pi)* 2*dreal( (8*CF*GS**2*NC*Pi*sqrtS**2*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*(muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    ((AR*muu1*(-1 + muu1**2 - muu2**2) + AL*muu2*(1 + muu1**2 - muu2**2))*dconjg(CL(i)) + (AL*muu1*(-1 + muu1**2 - muu2**2) + AR*muu2*(1 + muu1**2 - muu2**2))*dconjg(CR(i)))*
     &    ((BL*muu3*(-1 + muu3**2 - muu4**2) + BR*muu4*(1 + muu3**2 - muu4**2))*dconjg(DL(i)) + (BR*muu3*(-1 + muu3**2 - muu4**2) + BL*muu4*(1 + muu3**2 - muu4**2))*dconjg(DR(i))))/
     &  dreal(DenSmZ*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmH(i)))    )

** Z-Higgs dipole contribution 321

	  ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS -1/(4.*pi)* 2*dreal(  (8*CF*GS**2*NC*Pi*sqrtS**2*(-2*muu3**6*Sqrt(-4*muu3**2 + x1**2) + 2*muu4**6*Sqrt(-4*muu3**2 + x1**2) - 
     &      muu4**4*(-2*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) + 6*Sqrt(-4*muu3**2 + x1**2) + 
     &         x1*(Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) - 2*Sqrt(-4*muu3**2 + x1**2)) + Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2) + 
     &      muu3**4*(-6*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) + 3*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1 + 2*Sqrt(-4*muu3**2 + x1**2) - 
     &         2*muu4**2*Sqrt(-4*muu3**2 + x1**2) + 2*x1*Sqrt(-4*muu3**2 + x1**2) + 3*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2) + 
     &      muu4**2*(-(Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1**2) + 6*Sqrt(-4*muu3**2 + x1**2) + 
     &         2*x1*(Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) - 2*Sqrt(-4*muu3**2 + x1**2)) - 2*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2 + 
     &         Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2**2) + 
     &      (-1 + x1)*(2*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) + 2*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1**2 + 2*Sqrt(-4*muu3**2 + x1**2) - 
     &         3*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2 + Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2**2 + 
     &         Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1*(-5 + 3*x2)) + 
     &      muu3**2*(-8*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) + 14*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1 - 
     &         5*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1**2 + 2*Sqrt(-4*muu3**2 + x1**2) + 2*muu4**4*Sqrt(-4*muu3**2 + x1**2) - 4*x1*Sqrt(-4*muu3**2 + x1**2) + 
     &         6*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2 - 6*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1*x2 - 
     &         Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2**2 + 
     &         2*muu4**2*(-2*Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2)) + Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x1 - 2*Sqrt(-4*muu3**2 + x1**2) + 
     &            2*x1*Sqrt(-4*muu3**2 + x1**2) + Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*x2)))*
     &    ((AR*muu1*(-1 + muu1**2 - muu2**2) + AL*muu2*(1 + muu1**2 - muu2**2))*dconjg(CL(i)) + (AL*muu1*(-1 + muu1**2 - muu2**2) + AR*muu2*(1 + muu1**2 - muu2**2))*dconjg(CR(i)))*
     &    ((BL*muu3*(-1 + muu3**2 - muu4**2) + BR*muu4*(1 + muu3**2 - muu4**2))*dconjg(DL(i)) + (BR*muu3*(-1 + muu3**2 - muu4**2) + BL*muu4*(1 + muu3**2 - muu4**2))*dconjg(DR(i))))/
     &   dreal(DenSmZ*(-1 + muu3**2 + muu4**2)*(-1 - muu3**2 + muu4**2 + x1)**2*Sqrt(-4*muu3**2 + x1**2)*(-2 + x1 + x2)*dconjg(DenSmH(i)))   )

         endif


** Photon exchange in Chargino annihilation (m1=m2,m3=m4)
** Photon-Higgs interference

	if (icha.eq.2) then                                                     !Mo: H-G-Interferenz ausintegriert ist tatsächlich = 0! (Fuer 2->3 sowie fuer Dipol) (Grund??)
	
	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS + 0d0

        if (cutoff.lt.zerocheck) then

** G-Higgs dipole contribution 312

        ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - 0d0

** G-Higgs dipole contribution 321

        ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS -  0d0

        endif
		
	endif


       do j = jmin, jmax

       

** Higgs boson bremsstrahlung  33 + 44 + 34 + 43
        

        ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS +  (2*CF*gs**2*NC*sqrtS**2*(((-1 + muu1**2 + muu2**2)*CL(i) + 2*muu1*muu2*CR(i))*dconjg(CL(j)) 
     &     + (2*muu1*muu2*CL(i) + (-1 + muu1**2 + muu2**2)*CR(i))*dconjg(CR(j)))*
     &    (dconjg(DL(j))*((-2*muu3**6 + muu3**4*(2 + 2*muu4**2 + x1**2 - x2**2) + 
     &            muu3**2*(2 + 2*muu4**4 - x1**3 + (-2 + x2)**2*x2 + x1**2*(2 + x2) + x1*(-2 + x2)*(2 + 3*x2) - 4*muu4**2*(5 + x1**2 + 2*x1*(-2 + x2) + (-4 + x2)*x2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(-4 + x1**2 - (-2 + x2)*x2) - (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) + (-2 + x2)*x2)))*DL(i) - 
     &         4*muu3*muu4*(muu3**4 + muu3**2*(-2*(-1 + muu4**2 + x1) + (-2 + x1)*x2 + x2**2) + (1 + muu4**2 - x2)*(muu4**2 + (-1 + x1)*(-1 + x1 + x2)))*DR(i)) + 
     &      dconjg(DR(j))*(-4*muu3*muu4*(muu3**4 + muu3**2*(-2*(-1 + muu4**2 + x1) + (-2 + x1)*x2 + x2**2) + (1 + muu4**2 - x2)*(muu4**2 + (-1 + x1)*(-1 + x1 + x2)))*DL(i) + 
     &         (-2*muu3**6 + muu3**4*(2 + 2*muu4**2 + x1**2 - x2**2) + muu3**2*
     &             (2 + 2*muu4**4 - x1**3 + (-2 + x2)**2*x2 + x1**2*(2 + x2) + x1*(-2 + x2)*(2 + 3*x2) - 4*muu4**2*(5 + x1**2 + 2*x1*(-2 + x2) + (-4 + x2)*x2)) - 
     &            (1 + muu4**2 - x2)*(2*muu4**4 + muu4**2*(-4 + x1**2 - (-2 + x2)*x2) - (-1 + x1)*(2 + x1**2 + 2*x1*(-1 + x2) + (-2 + x2)*x2)))*DR(i))))/
     &  ((-1 - muu3**2 + muu4**2 + x1)**2*(-1 + muu3**2 - muu4**2 + x2)**2*dconjg(DenSmH(j))*DenSmH(i))

        if (cutoff.lt.zerocheck) then

** Hboson dipole contribution 312

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - (-2*CF*GS**2*NC*sqrtS**2*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*(muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (((-1 + muu1**2 + muu2**2)*CL(i) + 2*muu1*muu2*CR(i))*dconjg(CL(j)) + (2*muu1*muu2*CL(i) + (-1 + muu1**2 + muu2**2)*CR(i))*dconjg(CR(j)))*
     &    (dconjg(DL(j))*((-1 + muu3**2 + muu4**2)*DL(i) + 2*muu3*muu4*DR(i)) + dconjg(DR(j))*(2*muu3*muu4*DL(i) + (-1 + muu3**2 + muu4**2)*DR(i))))/
     &  ((-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmH(j))*DenSmH(i))

** Hboson dipole contribution 321

	ChiChi2QQ_RadiationS = ChiChi2QQ_RadiationS - (2*CF*GS**2*NC*sqrtS**2*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) + 
     &      (Sqrt((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*(-3*muu3**4 + muu4**4 + muu4**2*(x1 - x2) - (-1 + x1)*(-1 + 2*x1 + x2) + muu3**2*(-4 - 2*muu4**2 + 5*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (((-1 + muu1**2 + muu2**2)*CL(i) + 2*muu1*muu2*CR(i))*dconjg(CL(j)) + (2*muu1*muu2*CL(i) + (-1 + muu1**2 + muu2**2)*CR(i))*dconjg(CR(j)))*
     &    (dconjg(DL(j))*((-1 + muu3**2 + muu4**2)*DL(i) + 2*muu3*muu4*DR(i)) + dconjg(DR(j))*(2*muu3*muu4*DL(i) + (-1 + muu3**2 + muu4**2)*DR(i))))/
     &  ((1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmH(j))*DenSmH(i))
 
        endif

  
    

	  enddo


	enddo


	end

************************************************

	double precision function ChiChi2QQ_RadiationT()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j

        double precision temp1,temp2
	
***** itypeT,U follows the type of the sfermion exchanged in t,u-channel interactions
***** att gives the type of the antiquark final state
***** antitt is the opposite type to the final state quark type

	if (icha.eq.0) then
		itypeT = tt
		itypeU = tt
	else if (icha.eq.1) then
* tt is defined as normally as down-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt
		else if (tt.eq.4) then 
		      itypeT = tt
		      itypeU = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt+1
		else if (tt.eq.4) then 
		      itypeT = tt-1
		      itypeU = tt-1
		endif	
	endif	

	ChiChi2QQ_RadiationT = 0d0

	do j = 1, 2

	  do i = 1, 2

** Tchannel bremsstrahlung 11                                                       

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT + (-8*CF*gs**2*NC*(muu3**2*(dconjg(ELc(indi,j,tt,gen))*
     &          ((-2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*ELc(indi,i,tt,gen) + 2*muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)) + 
     &         dconjg(ER(indi,j,tt,gen))*(2*muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + 
     &            (-2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*ER(indi,i,tt,gen))) - 
     &      (-1 + muu3**2 - muu4**2 + x2)*(muu1*muu3*sqrtS**2*
     &          (dconjg(ER(indi,j,tt,gen))*ELc(indi,i,tt,gen) + dconjg(ELc(indi,j,tt,gen))*ER(indi,i,tt,gen)) + 
     &         ((2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*
     &            (dconjg(ELc(indi,j,tt,gen))*ELc(indi,i,tt,gen) + dconjg(ER(indi,j,tt,gen))*ER(indi,i,tt,gen)))/2.))*
     &    (muu2*muu4*sqrtS**2*(dconjg(FR(indj,j,tt,gen))*FL(indj,i,tt,gen) + dconjg(FL(indj,j,tt,gen))*FR(indj,i,tt,gen)) - 
     &      ((2*p1k2 - sqrtS**2*x2)*(dconjg(FL(indj,j,tt,gen))*FL(indj,i,tt,gen) + dconjg(FR(indj,j,tt,gen))*FR(indj,i,tt,gen)))/2.))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)**2*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))
   
** Tchannel bremsstrahlung 22

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT + (-8*CF*gs**2*NC*(dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen)) + 
     &      dconjg(ELc(indi,j,tt,gen))*(p1k1*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)))*
     &    (muu4**2*(dconjg(FL(indj,j,tt,gen))*((2*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1))*FL(indj,i,tt,gen) + 
     &            2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &         dconjg(FR(indj,j,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + 
     &            (2*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1))*FR(indj,i,tt,gen))) + 
     &      (1 + muu3**2 - muu4**2 - x1)*(muu2*muu4*sqrtS**2*
     &          (dconjg(FR(indj,j,tt,gen))*FL(indj,i,tt,gen) + dconjg(FL(indj,j,tt,gen))*FR(indj,i,tt,gen)) - 
     &         ((2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*
     &            (dconjg(FL(indj,j,tt,gen))*FL(indj,i,tt,gen) + dconjg(FR(indj,j,tt,gen))*FR(indj,i,tt,gen)))/2.)))/
     &  (sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)**2*(2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))

** Tchannel bremsstrahlung 33

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT +  (-4*CF*gs**2*NC*(2*p1k1 - 2*p1k2 - sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - x2))*
     &    (dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen)) + 
     &      dconjg(ELc(indi,j,tt,gen))*(p1k1*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)))*
     &    (dconjg(FL(indj,j,tt,gen))*((2*p1k2 - sqrtS**2*x2)*FL(indj,i,tt,gen) - 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) - 
     &      dconjg(FR(indj,j,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + (-2*p1k2 + sqrtS**2*x2)*FR(indj,i,tt,gen))))/
     &  ((2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*(2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*(2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))
     &    )
   
** Tchannel bremsstrahlung 12

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT + 2d0*dreal(  (CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(sqrtS**2*dconjg(FL(indj,j,tt,gen))*
     &          (2*muu1*muu3*ER(indi,i,tt,gen)*(sqrtS**2*
     &                (2*p1k2 - 2*muu3**2*p1k2 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu3**2*sqrtS**2 - 
     &                  muu1**2*muu3**2*sqrtS**2 + muu2**2*muu3**2*sqrtS**2 - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 - 
     &                  2*p1k1*x2 - 2*sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + muu3**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FL(indj,i,tt,gen) + 
     &               2*muu2*muu4*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*((2*p1k1*(-4*(-1 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-2 - muu1**2 + muu2**2 + muu3**2)*x2 + muu4**2*(2 + 2*muu1**2 - 2*muu2**2 + x2))) - 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 
     &                     (muu2**2 + muu4**2)*(-2*p1k2*x1 + (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2)) + 
     &                     muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu1**2*(2*p1k2*x1 + sqrtS**2*
     &                         (-1 + muu4**4 + x1 - 2*muu2**2*(1 + muu4**2 - x1 - x2) + muu4**2*(3*x1 - x2) + x2 - 2*x1*x2)) + 
     &                     muu3**2*(-2*p1k2*(-2 + 2*muu1**2 - 2*muu2**2 + x1) + 
     &                        sqrtS**2*(muu1**4 - muu1**2*(2*muu2**2 + 2*muu4**2 + x1 - x2) + 
     &                           (-1 + muu2**2)*(1 + muu2**2 - 2*muu4**2 - x1 + x2)))))*FL(indj,i,tt,gen) + 
     &               2*muu2*muu4*sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + 
     &                  muu4**2*sqrtS**2 + muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + 
     &                  muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*FR(indj,i,tt,gen))) + 
     &         dconjg(FR(indj,j,tt,gen))*(2*muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)*
     &             (2*muu2*muu4*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*FL(indj,i,tt,gen) + 
     &               sqrtS**2*(2*p1k2 - 2*muu3**2*p1k2 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu3**2*sqrtS**2 - 
     &                  muu1**2*muu3**2*sqrtS**2 + muu2**2*muu3**2*sqrtS**2 - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 - 
     &                  2*p1k1*x2 - 2*sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + muu3**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**4*
     &                (2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu4**2*sqrtS**2 + 
     &                  muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                  2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + 
     &                  muu2**2*sqrtS**2*x2)*FL(indj,i,tt,gen) + 
     &               sqrtS**2*(-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*
     &                   (muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu1**2*(1 + muu3**2 - x1 - x2) + 3*x2 - 2*x1*x2 - x2**2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2))))*FR(indj,i,tt,gen)))) + 
     &      dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     &          (2*muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen)*
     &             (sqrtS**2*(2*p1k2 - 2*muu3**2*p1k2 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu3**2*sqrtS**2 - 
     &                  muu1**2*muu3**2*sqrtS**2 + muu2**2*muu3**2*sqrtS**2 - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 - 
     &                  2*p1k1*x2 - 2*sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + muu3**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FL(indj,i,tt,gen) + 
     &               2*muu2*muu4*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*FR(indj,i,tt,gen)) + 
     &            ER(indi,i,tt,gen)*(sqrtS**2*(-4*p1k1**2*x2 + 
     &                  (2*p1k2 - sqrtS**2*x2)*(muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu1**2*(1 + muu3**2 - x1 - x2) + 3*x2 - 2*x1*x2 - x2**2 + 
     &                        muu2**2*(-1 - muu3**2 + x1 + x2))))*FL(indj,i,tt,gen) + 
     &               2*muu2*muu4*sqrtS**4*(2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + 
     &                  muu4**2*sqrtS**2 + muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + 
     &                  muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*FR(indj,i,tt,gen))) + 
     &         sqrtS**2*dconjg(FR(indj,j,tt,gen))*(-2*muu1*muu3*ELc(indi,i,tt,gen)*
     &             (2*muu2*muu4*sqrtS**4*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FL(indj,i,tt,gen) + 
     &               sqrtS**2*(-2*p1k2 + 2*muu3**2*p1k2 - sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 - muu3**2*sqrtS**2 + 
     &                  muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 
     &                  2*p1k1*x2 + 2*sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 - muu3**2*sqrtS**2*x2 - sqrtS**2*x1*x2 - 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FR(indj,i,tt,gen)) + 
     &            ER(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*
     &                (2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu4**2*sqrtS**2 + 
     &                  muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                  2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + 
     &                  muu2**2*sqrtS**2*x2)*FL(indj,i,tt,gen) + 
     &               (2*p1k1*(-4*(-1 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-2 - muu1**2 + muu2**2 + muu3**2)*x2 + muu4**2*(2 + 2*muu1**2 - 2*muu2**2 + x2))) - 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 
     &                     (muu2**2 + muu4**2)*(-2*p1k2*x1 + (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2)) + 
     &                     muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu1**2*(2*p1k2*x1 + sqrtS**2*
     &                         (-1 + muu4**4 + x1 - 2*muu2**2*(1 + muu4**2 - x1 - x2) + muu4**2*(3*x1 - x2) + x2 - 2*x1*x2)) + 
     &                     muu3**2*(-2*p1k2*(-2 + 2*muu1**2 - 2*muu2**2 + x1) + 
     &                        sqrtS**2*(muu1**4 - muu1**2*(2*muu2**2 + 2*muu4**2 + x1 - x2) + 
     &                           (-1 + muu2**2)*(1 + muu2**2 - 2*muu4**2 - x1 + x2)))))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))   )

** Tchannel bremsstrahlung 13

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT + 2d0*dreal( (2*CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*((muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + 
     &            p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + muu1**2*sqrtS**4*(1 + muu4**2 - x2))*
     &          ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + 
     &            sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*ER(indi,i,tt,gen)) + 
     &      dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*
     &          (-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*ELc(indi,i,tt,gen) + 
     &         (muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + 
     &            p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + muu1**2*sqrtS**4*(1 + muu4**2 - x2))*
     &          ER(indi,i,tt,gen)))*(dconjg(FL(indj,j,tt,gen))*
     &       ((-2*p1k2 + sqrtS**2*x2)*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &      dconjg(FR(indj,j,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + (-2*p1k2 + sqrtS**2*x2)*FR(indj,i,tt,gen))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))  )

** Tchannel bremsstrahlung 23

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT + 2d0*dreal(  (-2*CF*gs**2*NC*(dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen)) + 
     &      dconjg(ELc(indi,j,tt,gen))*(p1k1*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)))*
     &    (dconjg(FL(indj,j,tt,gen))*((4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + 
     &            sqrtS**2*(-2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + 
     &               sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + x1)*x2 + muu4**2*(2 - 2*muu1**2 - 2*x1 + x2))
     &               ))*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*
     &          (-2*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*FR(indj,i,tt,gen)) + 
     &      dconjg(FR(indj,j,tt,gen))*(-2*muu2*muu4*sqrtS**2*
     &          (2*p1k1 - 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*FL(indj,i,tt,gen) + 
     &         (4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + 
     &            sqrtS**2*(-2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + 
     &               sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + x1)*x2 + muu4**2*(2 - 2*muu1**2 - 2*x1 + x2))
     &               ))*FR(indj,i,tt,gen))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))   )

        if (cutoff.lt.zerocheck) then
             
** Tchannel dipole 312

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT -  (4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (dconjg(ELc(indi,j,tt,gen))*((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &            ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &            (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &            ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &             (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)) + 
     &      dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + 
     &         (((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &            ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &            (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &            ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &             (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*ER(indi,i,tt,gen)))*
     &    (dconjg(FL(indj,j,tt,gen))*(((muu1**2 - muu2**2 + muu3**2 - muu4**2)*sqrtS**2 - 
     &            2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &               (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &                (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FL(indj,i,tt,gen) - 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) - 
     &      dconjg(FR(indj,j,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + 
     &         ((-muu1**2 + muu2**2 - muu3**2 + muu4**2)*sqrtS**2 + 
     &            2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &               (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &                (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FR(indj,i,tt,gen))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*((muu1**2 + muu3**2)*sqrtS**2 - 
     &      2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &         ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &         (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &          (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    ((muu1**2 + muu3**2)*sqrtS**2 - 2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &         ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &         (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/
     &          (4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))

** Tchannel dipole 321

	ChiChi2QQ_RadiationT = ChiChi2QQ_RadiationT - (-4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (3*muu3**4 - muu4**4 + muu3**2*(4 + 2*muu4**2 - 5*x1 - x2) + muu4**2*(-x1 + x2) + (-1 + x1)*(-1 + 2*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (dconjg(ELc(indi,j,tt,gen))*((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &            ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &            (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &            ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &             (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)) + 
     &      dconjg(ER(indi,j,tt,gen))*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + 
     &         (((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &            ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &            (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &            ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &             (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*ER(indi,i,tt,gen)))*
     &    (dconjg(FL(indj,j,tt,gen))*(((muu1**2 - muu2**2 + muu3**2 - muu4**2)*sqrtS**2 - 
     &            2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &               (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &                (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FL(indj,i,tt,gen) - 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) - 
     &      dconjg(FR(indj,j,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + 
     &         ((-muu1**2 + muu2**2 - muu3**2 + muu4**2)*sqrtS**2 + 
     &            2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &               (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &                (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FR(indj,i,tt,gen))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*((muu1**2 + muu3**2)*sqrtS**2 - 
     &      2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &         ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &          (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    ((muu1**2 + muu3**2)*sqrtS**2 - 2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &         ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/
     &          (4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))

         endif


	  enddo

	enddo

	end
************************************************

	double precision function ChiChi2QQ_RadiationU()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j
        double precision temp1,temp2

***** itypeT,U follows the type of the sfermion exchanged in t,u-channel interactions
***** att gives the type of the antiquark final state
***** antitt is the opposite type to the final state quark type

	if (icha.eq.0) then
		itypeT = tt
		itypeU = tt
	else if (icha.eq.1) then
* tt is defined as normally as down-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt
		else if (tt.eq.4) then 
		      itypeT = tt
		      itypeU = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt+1
		else if (tt.eq.4) then 
		      itypeT = tt-1
		      itypeU = tt-1
		endif	
	endif	

	ChiChi2QQ_RadiationU = 0d0

	do j = 1, 2

	  do i = 1, 2

** Uchannel bremsstrahlung 11                                                        

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU +  (-8*CF*gs**2*NC*(dconjg(GR(indi,j,tt,gen))*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen)) + 
     &      dconjg(GL(indi,j,tt,gen))*(p1k2*GL(indi,i,tt,gen) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)))*
     &    (muu3**2*(dconjg(HL(indj,j,tt,gen))*((2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x2))*HL(indj,i,tt,gen) + 
     &            2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &         dconjg(HR(indj,j,tt,gen))*(2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + 
     &            (2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x2))*HR(indj,i,tt,gen))) - 
     &      (-1 + muu3**2 - muu4**2 + x2)*(muu2*muu3*sqrtS**2*
     &          (dconjg(HR(indj,j,tt,gen))*HL(indj,i,tt,gen) + dconjg(HL(indj,j,tt,gen))*HR(indj,i,tt,gen)) - 
     &         ((2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*
     &            (dconjg(HL(indj,j,tt,gen))*HL(indj,i,tt,gen) + dconjg(HR(indj,j,tt,gen))*HR(indj,i,tt,gen)))/2.)))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)**2*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))

** Uchannel bremsstrahlung 22

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU + (-8*CF*gs**2*NC*(muu4**2*(dconjg(GL(indi,j,tt,gen))*
     &          ((-2*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2)*GL(indi,i,tt,gen) + 2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)) + 
     &         dconjg(GR(indi,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + 
     &            (-2*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2)*GR(indi,i,tt,gen))) + 
     &      (1 + muu3**2 - muu4**2 - x1)*(muu1*muu4*sqrtS**2*
     &          (dconjg(GR(indi,j,tt,gen))*GL(indi,i,tt,gen) + dconjg(GL(indi,j,tt,gen))*GR(indi,i,tt,gen)) + 
     &         ((2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*
     &            (dconjg(GL(indi,j,tt,gen))*GL(indi,i,tt,gen) + dconjg(GR(indi,j,tt,gen))*GR(indi,i,tt,gen)))/2.))*
     &    (muu2*muu3*sqrtS**2*(dconjg(HR(indj,j,tt,gen))*HL(indj,i,tt,gen) + dconjg(HL(indj,j,tt,gen))*HR(indj,i,tt,gen)) - 
     &      ((2*p1k1 - sqrtS**2*x1)*(dconjg(HL(indj,j,tt,gen))*HL(indj,i,tt,gen) + dconjg(HR(indj,j,tt,gen))*HR(indj,i,tt,gen)))/2.))/
     &  (sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)**2*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))

** Uchannel bremsstrahlung 33

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU + (4*CF*gs**2*NC*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu1**2 + muu2**2 + muu3**2 + muu4**2 - x1))*
     &    (dconjg(GR(indi,j,tt,gen))*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen)) + 
     &      dconjg(GL(indi,j,tt,gen))*(p1k2*GL(indi,i,tt,gen) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)))*
     &    (dconjg(HL(indj,j,tt,gen))*((2*p1k1 - sqrtS**2*x1)*HL(indj,i,tt,gen) - 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) - 
     &      dconjg(HR(indj,j,tt,gen))*(2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + (-2*p1k1 + sqrtS**2*x1)*HR(indj,i,tt,gen))))/
     &  ((2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen)))*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen)))
     &    )

** Uchannel bremsstrahlung 12

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU + 2*dreal(  (CF*gs**2*NC*(dconjg(GL(indi,j,tt,gen))*(sqrtS**2*dconjg(HL(indj,j,tt,gen))*
     &          (2*muu1*muu4*GR(indi,i,tt,gen)*(sqrtS**2*
     &                (-2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu4**2*sqrtS**2 - 
     &                  muu1**2*muu4**2*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 2*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - 
     &                  muu2**2*sqrtS**2*x1 + muu4**2*sqrtS**2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) - 
     &                  sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2)*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((2*p1k1*(-4*(-1 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-muu1**2 + muu2**2 + muu3**2)*x2 + muu4**2*(-2 + 2*muu1**2 - 2*muu2**2 + x2))) - 
     &                  sqrtS**2*(-(muu2**2*sqrtS**2) + muu2**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 - muu4**2*sqrtS**2 + 
     &                     muu2**4*muu4**2*sqrtS**2 - muu4**4*sqrtS**2 + muu2**2*muu4**4*sqrtS**2 + 4*p1k2*x1 - 2*muu2**2*p1k2*x1 - 
     &                     2*muu4**2*p1k2*x1 + muu2**2*sqrtS**2*x1 - muu2**4*sqrtS**2*x1 - muu4**2*sqrtS**2*x1 + 
     &                     muu2**2*muu4**2*sqrtS**2*x1 - 
     &                     muu3**2*(2*p1k2*(2 + 2*muu1**2 - 2*muu2**2 + x1) - 
     &                        sqrtS**2*(muu1**4 - muu1**2*(2*muu2**2 + 2*muu4**2 + x1 - 3*x2) + 
     &                           (-1 + muu2**2)*(1 + muu2**2 - 2*muu4**2 - x1 - x2))) + muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu2**2*sqrtS**2*x2 - muu2**4*sqrtS**2*x2 + muu4**2*sqrtS**2*x2 - muu2**2*muu4**2*sqrtS**2*x2 + 
     &                     muu1**2*(2*p1k2*x1 + sqrtS**2*
     &                         (-1 + muu4**4 + x1 - 2*muu2**2*(1 + muu4**2 - x1 - x2) + muu4**2*(x1 - x2) + x2 - 2*x1*x2))))*
     &                HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*
     &                (-2*p1k2 + 2*muu3**2*p1k2 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + 
     &                  muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                  muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HR(indj,i,tt,gen))) + 
     &         dconjg(HR(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*
     &             (2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HL(indj,i,tt,gen) + 
     &               sqrtS**2*(-2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu4**2*sqrtS**2 - 
     &                  muu1**2*muu4**2*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 2*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - 
     &                  muu2**2*sqrtS**2*x1 + muu4**2*sqrtS**2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) - 
     &                  sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**4*
     &                (-2*p1k2 + 2*muu3**2*p1k2 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + 
     &                  muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                  muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HL(indj,i,tt,gen) + 
     &               sqrtS**2*(8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 
     &                  2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 
     &                  6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 
     &                  2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + 
     &                  muu4**2*(8*p1k1**2 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1) + 
     &                     sqrtS**2*((2 - 2*muu1**2 + 2*muu2**2)*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)) - 4*p1k1**2*x2 - 
     &                  2*p1k2*sqrtS**2*x2 + 2*muu1**2*p1k2*sqrtS**2*x2 - 2*muu2**2*p1k2*sqrtS**2*x2 - 4*muu3**2*p1k2*sqrtS**2*x2 + 
     &                  4*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 - muu1**2*sqrtS**4*x1*x2 + muu2**2*sqrtS**4*x1*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*
     &                      (1 + muu3**2 - x1 + muu1**2*(1 + muu3**2 - x1 - x2) - x2 - x1*x2 + muu2**2*(-1 - muu3**2 + x1 + x2))))*
     &                HR(indj,i,tt,gen)))) + dconjg(GR(indi,j,tt,gen))*
     &       (dconjg(HL(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             (sqrtS**2*(-2*(-1 + muu3**2 + muu4**2)*p1k1 + sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu4**2*sqrtS**2 - 
     &                  muu1**2*muu4**2*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 2*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - 
     &                  muu2**2*sqrtS**2*x1 + muu4**2*sqrtS**2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) - 
     &                  sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2)*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(sqrtS**2*(8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 
     &                  2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 
     &                  6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 
     &                  2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + 
     &                  muu4**2*(8*p1k1**2 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1) + 
     &                     sqrtS**2*((2 - 2*muu1**2 + 2*muu2**2)*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)) - 4*p1k1**2*x2 - 
     &                  2*p1k2*sqrtS**2*x2 + 2*muu1**2*p1k2*sqrtS**2*x2 - 2*muu2**2*p1k2*sqrtS**2*x2 - 4*muu3**2*p1k2*sqrtS**2*x2 + 
     &                  4*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 - muu1**2*sqrtS**4*x1*x2 + muu2**2*sqrtS**4*x1*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*
     &                      (1 + muu3**2 - x1 + muu1**2*(1 + muu3**2 - x1 - x2) - x2 - x1*x2 + muu2**2*(-1 - muu3**2 + x1 + x2))))*
     &                HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**4*
     &                (-2*p1k2 + 2*muu3**2*p1k2 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + 
     &                  muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                  muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HR(indj,i,tt,gen))) + 
     &         sqrtS**2*dconjg(HR(indj,j,tt,gen))*(-2*muu1*muu4*GL(indi,i,tt,gen)*
     &             (2*muu2*muu3*sqrtS**4*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HL(indj,i,tt,gen) + 
     &               sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 - sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 - muu4**2*sqrtS**2 + 
     &                  muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + 2*p1k2*x1 + 2*sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 - muu4**2*sqrtS**2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) + 
     &                  sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 - sqrtS**2*x1*x2)*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*
     &                (-2*p1k2 + 2*muu3**2*p1k2 + sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + 
     &                  muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                  muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + 
     &                  muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HL(indj,i,tt,gen) + 
     &               (2*p1k1*(-4*(-1 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-muu1**2 + muu2**2 + muu3**2)*x2 + muu4**2*(-2 + 2*muu1**2 - 2*muu2**2 + x2))) - 
     &                  sqrtS**2*(-(muu2**2*sqrtS**2) + muu2**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 - muu4**2*sqrtS**2 + 
     &                     muu2**4*muu4**2*sqrtS**2 - muu4**4*sqrtS**2 + muu2**2*muu4**4*sqrtS**2 + 4*p1k2*x1 - 2*muu2**2*p1k2*x1 - 
     &                     2*muu4**2*p1k2*x1 + muu2**2*sqrtS**2*x1 - muu2**4*sqrtS**2*x1 - muu4**2*sqrtS**2*x1 + 
     &                     muu2**2*muu4**2*sqrtS**2*x1 - 
     &                     muu3**2*(2*p1k2*(2 + 2*muu1**2 - 2*muu2**2 + x1) - 
     &                        sqrtS**2*(muu1**4 - muu1**2*(2*muu2**2 + 2*muu4**2 + x1 - 3*x2) + 
     &                           (-1 + muu2**2)*(1 + muu2**2 - 2*muu4**2 - x1 - x2))) + muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu2**2*sqrtS**2*x2 - muu2**4*sqrtS**2*x2 + muu4**2*sqrtS**2*x2 - muu2**2*muu4**2*sqrtS**2*x2 + 
     &                     muu1**2*(2*p1k2*x1 + sqrtS**2*
     &                         (-1 + muu4**4 + x1 - 2*muu2**2*(1 + muu4**2 - x1 - x2) + muu4**2*(x1 - x2) + x2 - 2*x1*x2))))*
     &                HR(indj,i,tt,gen))))))/
     &  (sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*
     &    (-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))  )

** Uchannel bremsstrahlung 13

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU + 2*dreal( (2*CF*gs**2*NC*(dconjg(GR(indi,j,tt,gen))*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen)) + 
     &      dconjg(GL(indi,j,tt,gen))*(p1k2*GL(indi,i,tt,gen) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)))*
     &    (dconjg(HL(indj,j,tt,gen))*((-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &            sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + 
     &               sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*HL(indj,i,tt,gen) + 
     &         2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*HR(indj,i,tt,gen)
     &         ) + dconjg(HR(indj,j,tt,gen))*(-2*muu2*muu3*sqrtS**2*
     &          (-2*p1k1 + 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*HL(indj,i,tt,gen) + 
     &         (-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &            sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + 
     &               sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*HR(indj,i,tt,gen))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen)))) )

** Uchannel bremsstrahlung 23

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU + 2*dreal(  (-2*CF*gs**2*NC*(dconjg(GR(indi,j,tt,gen))*(muu1*muu4*sqrtS**2*
     &          (2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*GL(indi,i,tt,gen) + 
     &         (p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + sqrtS**2*
     &             (sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + 
     &               p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*GR(indi,i,tt,gen)) + 
     &      dconjg(GL(indi,j,tt,gen))*((p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + 
     &            sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + 
     &               p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*GL(indi,i,tt,gen) + 
     &         muu1*muu4*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*GR(indi,i,tt,gen)))*
     &    (dconjg(HL(indj,j,tt,gen))*((-2*p1k1 + sqrtS**2*x1)*HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &      dconjg(HR(indj,j,tt,gen))*(2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + (-2*p1k1 + sqrtS**2*x1)*HR(indj,i,tt,gen))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))  )

        if (cutoff.lt.zerocheck) then


** Uchannel dipole 312

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU -  (-4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (dconjg(GL(indi,j,tt,gen))*(((muu1**2 - muu2**2 - muu3**2 + muu4**2)*sqrtS**2 + 
     &            (-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - 
     &                  (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                  x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                     (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/(2.*sqrtS**4*(4*muu4**2 - x2**2)))*GL(indi,i,tt,gen) + 
     &         2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)) + 
     &      dconjg(GR(indi,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + 
     &         ((muu1**2 - muu2**2 - muu3**2 + muu4**2)*sqrtS**2 + 
     &            (-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - 
     &                  (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                  x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                     (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/(2.*sqrtS**4*(4*muu4**2 - x2**2)))*GR(indi,i,tt,gen)))*
     &    (dconjg(HL(indj,j,tt,gen))*(((-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &               Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &              sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - 
     &                 (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                 x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                    (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                     Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HL(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu4**2 - x2**2)) + 
     &         muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + dconjg(HR(indj,j,tt,gen))*
     &       (muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + ((-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &               Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &              sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - 
     &                 (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                 x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                    (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                     Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HR(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu4**2 - x2**2)))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         sqrtS**2*(4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**4 + (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &            x2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2) + 
     &               (-1 - muu1**2 + muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/(2.*sqrtS**4*(4*muu4**2 - x2**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    ((muu2**2 + muu3**2)*sqrtS**2 + (4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &          Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + sqrtS**2*
     &          (4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**4 + (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &            x2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2) + 
     &               (-1 - muu1**2 + muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/(2.*sqrtS**4*(4*muu4**2 - x2**2)) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))

** Uchannel dipole 321

	ChiChi2QQ_RadiationU = ChiChi2QQ_RadiationU -  (4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     &      (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (3*muu3**4 - muu4**4 + muu3**2*(4 + 2*muu4**2 - 5*x1 - x2) + muu4**2*(-x1 + x2) + (-1 + x1)*(-1 + 2*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (dconjg(GL(indi,j,tt,gen))*(-((4*(1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + 
     &               (1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &               4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &               (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1*
     &                Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - (1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*
     &             GL(indi,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) + 2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)) + 
     &      dconjg(GR(indi,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) - 
     &         ((4*(1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &              4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1*
     &               Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - (1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*
     &            GR(indi,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2))))*
     &    (dconjg(HL(indj,j,tt,gen))*(-((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &               sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - 
     &                  (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                  x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &                     (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*HL(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu3**2 - x1**2)) + 
     &         muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + dconjg(HR(indj,j,tt,gen))*
     &       (muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) - ((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &               Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &              sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - 
     &                 (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                 x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &                    (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                     Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*HR(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu3**2 - x1**2)))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &         sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - 
     &            (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &            x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &               (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    ((muu2**2 + muu3**2)*sqrtS**2 + (-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &          Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + sqrtS**2*
     &          (4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &            x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &               (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*
     &                Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) - (MSf2(j,itypeU,gen)+cI*MSf(j,itypeU,gen)*WSf(j,itypeU,gen))))

          endif


	  enddo

	enddo

	end


************************************************

	double precision function ChiChi2QQ_RadiationTS()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j,imin,imax,jmin,jmax
	complex*16 DenSmZ, DenSmH(6)
	double precision DenSmZ1, DenSmH1(6)

c       wZ = 0d0
c	wW = 0d0
	
c	wHk0(1) = 0d0
c	wHk0(2) = 0d0
c	wHk0(3) = 0d0
c	wHk0(4) = 0d0

c	wHk(1) = 0d0
c	wHk(2) = 0d0
c	wHk(3) = 0d0
c	wHk(4) = 0d0
c	wHk(5) = 0d0
c	wHk(6) = 0d0

	if (icha.eq.0) then
		itypeT = tt
		itypeU = tt
	else if (icha.eq.1) then
* tt is defined as normally as down-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt
		else if (tt.eq.4) then 
		      itypeT = tt
		      itypeU = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt+1
		else if (tt.eq.4) then 
		      itypeT = tt-1
		      itypeU = tt-1
		endif	
	endif		
    

        if (icha.eq.1) then
c only charged Higgs exchange		
	  imin = 5
	  imax = 6 
	  jmin = 5
	  jmax = 6
	else
c exchange of neutral Higgses 		
	  imin = 1
	  imax = 4 
	  jmin = 1
	  jmax = 4	  	
	endif

	if (icha.eq.1) then

	DenSmZ = (sqrtS**2 - MW2) + cI*MW*wW
        DenSmZc = (sqrtS**2 - MW2) - cI*MW*wW
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MW2)
		
	else
		
	DenSmZ = (sqrtS**2 - MZ2) + cI*MZ*wZ
        DenSmZc = (sqrtS**2 - MZ2) - cI*MZ*wZ
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MZ2)

	endif

	do i = imin, imax
	  DenSmH(i) = (sqrtS**2 - MHk0(i)**2) + cI*MHk0(i)*wHk0(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk0(i)**2)
	enddo

	ChiChi2QQ_RadiationTS = 0d0

	do i = 1, 2


** TSchannel bremsstrahlung 11

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(-4*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)*((4*p1k2**2*(-1 + 2*muu3**2 + x2) - 2*p1k2*sqrtS**2*(1 + x2)*(-1 + 2*muu3**2 + x2) 
     &           - 2*p1k1*(-1 + x2)*(-2*p1k2 + sqrtS**2*x2) + 
     &            sqrtS**4*(2*muu3**4 + (-1 + x2)*x2 + muu3**2*(1 + x2) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*FL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(-2*p1k1*(-1 + x2) - sqrtS**2*(-1 + x2)*(-3 + x1 + 3*x2) + muu3**2*(-8*p1k2 + 2*sqrtS**2*(3 + x2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(-2*p1k1*(-1 + x2) - sqrtS**2*(-1 + x2)*(-3 + x1 + 3*x2) + muu3**2*(-8*p1k2 + 2*sqrtS**2*(3 + x2)))*FL(indj,i,tt,gen) + 
     &         (4*p1k2**2*(-1 + 2*muu3**2 + x2) - 2*p1k2*sqrtS**2*(1 + x2)*(-1 + 2*muu3**2 + x2) - 2*p1k1*(-1 + x2)*(-2*p1k2 + sqrtS**2*x2) + 
     &            sqrtS**4*(2*muu3**4 + (-1 + x2)*x2 + muu3**2*(1 + x2) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*FR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x2)**2*(-2*p1k2 - sqrtS**2*(muu1**2 + muu3**2 - x2) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

	endif	


	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(   -(4*CF*gs**2*NC*(-(dconjg(AR)*(-(dconjg(BL)*(2*muu1*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*ELc(indi,i,tt,gen) - 
     &                (2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*ER(indi,i,tt,gen))*
     &              ((-2*p1k2 + sqrtS**2*x2)*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen))) + 
     &           sqrtS**2*dconjg(BR)*(muu3*(1 + muu3**2 + muu4**2 - x2)*ER(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FL(indj,i,tt,gen) - 2*muu2*p1k2*FR(indj,i,tt,gen)) + 
     &              muu1*ELc(indi,i,tt,gen)*(muu4*((1 - muu1**2 + muu2**2)*sqrtS**2 + 2*p1k2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*FL(indj,i,tt,gen) - 
     &                 muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*FR(indj,i,tt,gen))))) + 
     &      dconjg(AL)*(dconjg(BR)*(-((2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*ELc(indi,i,tt,gen)) + 
     &            2*muu1*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*ER(indi,i,tt,gen))*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + (-2*p1k2 + sqrtS**2*x2)*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu3*(1 + muu3**2 + muu4**2 - x2)*ELc(indi,i,tt,gen)*(2*muu2*p1k2*FL(indj,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*FL(indj,i,tt,gen) + 
     &               muu4*((-1 + muu1**2 - muu2**2)*sqrtS**2 + 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 2*p1k2*(-1 + 3*muu3**2 - muu4**2 + x2) - 
     &                  sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)**2*dconjg(DenSmZ)*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )


** TSchannel bremsstrahlung 12

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal( -(-2*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)*((-4*p1k2**2*x1 - 4*p1k1**2*x2 + 2*p1k2*sqrtS**2*(-1 + x2 + x1*(1 + 2*muu3**2 + x2)) + 
     &            sqrtS**4*(4*muu3**4 + x2 - x2*(x1 + x2) - 2*muu3**2*(-1 + 2*x1 + x2) + muu1**2*(4*muu3**2 - 2*(-1 + x1 + x2))) + 
     &            2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(1 + (-3 + x2)*x2 + 2*muu3**2*(2 + x2) + x1*(-1 + 2*x2))))*FL(indj,i,tt,gen) + 
     &         2*muu1*muu3*sqrtS**2*(2*p1k2*x1 + 2*p1k1*(-1 + 4*muu3**2 + x1 - x2) + sqrtS**2*(4 - 2*muu3**2*(-3 + x1) - 5*x1 + 2*(-2 + x1)*x2))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(2*muu1*muu3*sqrtS**2*(2*p1k2*x1 + 2*p1k1*(-1 + 4*muu3**2 + x1 - x2) + sqrtS**2*(4 - 2*muu3**2*(-3 + x1) - 5*x1 + 2*(-2 + x1)*x2))*FL(indj,i,tt,gen) + 
     &         (-4*p1k2**2*x1 - 4*p1k1**2*x2 + 2*p1k2*sqrtS**2*(-1 + x2 + x1*(1 + 2*muu3**2 + x2)) + 
     &            sqrtS**4*(4*muu3**4 + x2 - x2*(x1 + x2) - 2*muu3**2*(-1 + 2*x1 + x2) + muu1**2*(4*muu3**2 - 2*(-1 + x1 + x2))) + 
     &            2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(1 + (-3 + x2)*x2 + 2*muu3**2*(2 + x2) + x1*(-1 + 2*x2))))*FR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)*(-1 + x2)*(-2*p1k2 - sqrtS**2*(muu1**2 + muu3**2 - x2) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))  )

	endif	

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(2*CF*gs**2*NC*(dconjg(AR)*(dconjg(BL)*(-(ER(indi,i,tt,gen)*((-4*p1k1**2*x2 + 
     &                    (2*p1k2 - sqrtS**2*x2)*(muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                    muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                    2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2 - x2**2)))*FL(indj,i,tt,gen) + 
     &                 2*muu2*muu4*sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*FR(indj,i,tt,gen))) + 
     &            2*muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen)*(-((-2*(-1 + muu3**2)*p1k2 - 2*p1k1*x2 + 
     &                    sqrtS**2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1)) + (-2 + muu1**2 - muu2**2 + muu3**2 + x1)*x2) + 
     &                    muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FL(indj,i,tt,gen)) + 
     &               2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FR(indj,i,tt,gen))) + 
     &         dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(muu4*sqrtS**2*(2*p1k2*x1 + 2*p1k1*(-2 + x2) + 
     &                  sqrtS**2*(-1 - muu4**2 + 3*x1 + muu1**2*(1 + muu4**2 - x1 - x2) + x2 - 2*x1*x2 + muu2**2*(-1 - muu4**2 + x1 + x2)) + 
     &                  muu3**2*(-8*p1k2 + sqrtS**2*(-3 + 3*muu1**2 - 3*muu2**2 + 4*x2)))*FL(indj,i,tt,gen) + 
     &               2*muu2*sqrtS**4*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*(2*muu4*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 
     &                  2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(-2*p1k2*(-2 + 4*muu4**2 + x1) - (1 + muu1**2 - muu2**2)*sqrtS**2*(3*muu3**2 - 3*(-1 + muu4**2 + x1) - x2) + 
     &                  p1k1*(4*muu3**2 - 2*(2*(-1 + muu4**2 + x1) + x2)))*FR(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(dconjg(BL)*(muu3*ELc(indi,i,tt,gen)*(muu2*sqrtS**2*(-2*p1k2*(-2 + 4*muu4**2 + x1) - (1 + muu1**2 - muu2**2)*sqrtS**2*(3*muu3**2 - 3*(-1 + muu4**2 + x1) - x2) + 
     &                  p1k1*(4*muu3**2 - 2*(2*(-1 + muu4**2 + x1) + x2)))*FL(indj,i,tt,gen) + 
     &               2*muu4*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*FR(indj,i,tt,gen)) + 
     &            muu1*sqrtS**2*ER(indi,i,tt,gen)*(2*muu2*sqrtS**2*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*FL(indj,i,tt,gen) + 
     &               muu4*(2*p1k2*x1 + 2*p1k1*(-2 + x2) + sqrtS**2*(-1 - muu4**2 + 3*x1 + muu1**2*(1 + muu4**2 - x1 - x2) + x2 - 2*x1*x2 + muu2**2*(-1 - muu4**2 + x1 + x2)) + 
     &                  muu3**2*(-8*p1k2 + sqrtS**2*(-3 + 3*muu1**2 - 3*muu2**2 + 4*x2)))*FR(indj,i,tt,gen))) + 
     &         dconjg(BR)*(2*muu1*muu3*sqrtS**2*ER(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FL(indj,i,tt,gen) - 
     &               (-2*(-1 + muu3**2)*p1k2 - 2*p1k1*x2 + sqrtS**2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1)) + (-2 + muu1**2 - muu2**2 + muu3**2 + x1)*x2) + 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*(-2*muu2*muu4*sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*FL(indj,i,tt,gen) - 
     &               (-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*(muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2 - x2**2)))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmZ)*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))  )


** TSchannel bremsstrahlung 21

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(   -(2*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)*((-4*p1k2**2*x1 - 4*p1k1**2*x2 + sqrtS**4
     &            *(1 + 2*muu3**2 - x1 - x2)*(2*(muu1**2 + muu3**2) + x2) + 
     &            2*p1k2*sqrtS**2*(-1 + 2*muu3**2*(-2 + x1) + x1 + x2 + x1*x2) + 2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(1 + x2*(-3 + 2*muu3**2 + x2) + x1*(-1 + 2*x2))))*
     &          FL(indj,i,tt,gen) + 2*muu1*muu3*sqrtS**2*(2*p1k2*(1 - 4*muu3**2 + x1 - x2) - 2*p1k1*x2 + sqrtS**2*(4 + 2*x1*(-2 + x2) + (-6 + x2)*x2 + 2*muu3**2*(3 + x2)))*
     &          FR(indj,i,tt,gen)) + ELc(indi,i,tt,gen)*(2*muu1*muu3*sqrtS**2*
     &          (2*p1k2*(1 - 4*muu3**2 + x1 - x2) - 2*p1k1*x2 + sqrtS**2*(4 + 2*x1*(-2 + x2) + (-6 + x2)*x2 + 2*muu3**2*(3 + x2)))*FL(indj,i,tt,gen) + 
     &         (-4*p1k2**2*x1 - 4*p1k1**2*x2 + sqrtS**4*(1 + 2*muu3**2 - x1 - x2)*(2*(muu1**2 + muu3**2) + x2) + 2*p1k2*sqrtS**2*(-1 + 2*muu3**2*(-2 + x1) + x1 + x2 + x1*x2) + 
     &            2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(1 + x2*(-3 + 2*muu3**2 + x2) + x1*(-1 + 2*x2))))*FR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)*(-1 + x2)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

	endif	

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(2*CF*gs**2*NC*(sqrtS**2*dconjg(AR)*(dconjg(BL)*(-(ER(indi,i,tt,gen)*
     &               ((-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*(muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                    muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                    2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2 - x2**2)))*FL(indj,i,tt,gen) + 
     &                 2*muu2*muu4*sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*FR(indj,i,tt,gen))) + 
     &            2*muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen)*(-((-2*(-1 + muu3**2)*p1k2 - 2*p1k1*x2 + 
     &                    sqrtS**2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1)) + (-2 + muu1**2 - muu2**2 + muu3**2 + x1)*x2) + 
     &                    muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FL(indj,i,tt,gen)) + 
     &               2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FR(indj,i,tt,gen))) + 
     &         dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(muu4*sqrtS**2*(2*p1k1*(-2 + 4*muu3**2 + x2) + 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) - 
     &                  sqrtS**2*(3 + 3*muu4**2 - 3*x1 + muu3**2*(-3 + 4*x1) - 5*x2 + 2*(muu3**2 - muu4**2 + x1)*x2 + 2*x2**2 + muu1**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2) - 
     &                     muu2**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2)))*FL(indj,i,tt,gen) + 
     &               2*muu2*sqrtS**4*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*(2*muu4*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*
     &                   (muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                  2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu4**2*(8*p1k1 - 3*(1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*(-2 + x1) - 2*p1k1*x2 + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + x1 + x2))*
     &                FR(indj,i,tt,gen)))) + dconjg(AL)*(sqrtS**2*dconjg(BL)*
     &          (muu3*ELc(indi,i,tt,gen)*(muu2*sqrtS**2*(muu4**2*(8*p1k1 - 3*(1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*(-2 + x1) - 2*p1k1*x2 + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + x1 + x2))*FL(indj,i,tt,gen) + 
     &               2*muu4*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*(muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                  2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*FR(indj,i,tt,gen)) - 
     &            muu1*ER(indi,i,tt,gen)*(2*muu2*sqrtS**4*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))*FL(indj,i,tt,gen) + 
     &               muu4*sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                  sqrtS**2*(3 + 3*muu4**2 - 3*x1 + muu3**2*(-3 + 4*x1) - 5*x2 + 2*(muu3**2 - muu4**2 + x1)*x2 + 2*x2**2 + muu1**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2) - 
     &                     muu2**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2)))*FR(indj,i,tt,gen))) + 
     &         dconjg(BR)*(2*muu1*muu3*sqrtS**4*ER(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FL(indj,i,tt,gen) - 
     &               (-2*(-1 + muu3**2)*p1k2 - 2*p1k1*x2 + sqrtS**2*(-((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1)) + (-2 + muu1**2 - muu2**2 + muu3**2 + x1)*x2) + 
     &                  muu4**2*(4*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 + x2)))*FR(indj,i,tt,gen)) + 
     &            sqrtS**2*ELc(indi,i,tt,gen)*(-2*muu2*muu4*sqrtS**2*(2*(-1 + muu3**2 + muu4**2)*p1k1 + muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*FL(indj,i,tt,gen) - 
     &               (-4*p1k1**2*x2 + (2*p1k2 - sqrtS**2*x2)*(muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2)) + 
     &                  muu4**2*(8*p1k1**2 + 2*p1k1*sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)) - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2 - x2**2)))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmZ)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))) )


** TSchannel bremsstrahlung 22

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(4*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)*((4*p1k1**2*(-1 + 2*muu3**2 + x1) - 2*p1k1*(-1 + x1)
     &            *(-2*p1k2 + sqrtS**2*(-1 + 2*muu3**2 + x1 + x2)) + 
     &            sqrtS**4*(muu3**2*(1 + 2*muu3**2 - x1) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*FL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(muu3**2*(8*p1k1 - 2*sqrtS**2*(-3 + x1)) + (-1 + x1)*(2*p1k2 + sqrtS**2*(3 - 3*x1 - 2*x2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(muu3**2*(8*p1k1 - 2*sqrtS**2*(-3 + x1)) + (-1 + x1)*(2*p1k2 + sqrtS**2*(3 - 3*x1 - 2*x2)))*FL(indj,i,tt,gen) + 
     &         (4*p1k1**2*(-1 + 2*muu3**2 + x1) - 2*p1k1*(-1 + x1)*(-2*p1k2 + sqrtS**2*(-1 + 2*muu3**2 + x1 + x2)) + 
     &            sqrtS**4*(muu3**2*(1 + 2*muu3**2 - x1) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*FR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)**2*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

	endif	

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal( -(4*CF*gs**2*NC*(dconjg(AL)*(-2*dconjg(BR)*(p1k1*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen))*
     &          (-2*muu2*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FL(indj,i,tt,gen) + 2*p1k1*FR(indj,i,tt,gen) + 
     &            (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 + x1*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-2 + muu1**2 - muu2**2 + x1)) + sqrtS**2*(-1 + x1)*x2 + 
     &               muu3**2*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2)) + muu4**2*(-6*p1k1 - 2*p1k2 + sqrtS**2*(3*(-1 + muu1**2 - muu2**2 + x1) + x2)))*FR(indj,i,tt,gen)
     &            ) + sqrtS**2*dconjg(BL)*(muu3*ELc(indi,i,tt,gen)*(muu2*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + 3*muu4**2 + x1))*FL(indj,i,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*FL(indj,i,tt,gen) + 
     &               muu4*(1 + muu3**2 + muu4**2 - x1)*(-2*p1k1 + sqrtS**2*x1)*FR(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(2*dconjg(BL)*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen))*
     &          (-((2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                 sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*FL(indj,i,tt,gen)) + 
     &            2*muu2*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu3*ER(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FL(indj,i,tt,gen)) + 
     &               muu2*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + 3*muu4**2 + x1))*FR(indj,i,tt,gen))
     &             + muu1*ELc(indi,i,tt,gen)*(muu4*(1 + muu3**2 + muu4**2 - x1)*(-2*p1k1 + sqrtS**2*x1)*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)**2*dconjg(DenSmZ)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

** TSchannel bremsstrahlung 31

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal( -(-2*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)
     &            *((-4*muu3**4*sqrtS**6 + 2*p1k1*(8*p1k2**2 + 2*p1k2*sqrtS**2*(-1 + 2*muu3**2 - 3*x2) + sqrtS**4*x2*(1 + muu1**2 - muu3**2 + x2)) + 
     &            muu1**2*sqrtS**4*(-2*p1k2*(-2 + x1 + 2*x2) + sqrtS**2*x2*(-3 + x1 + 3*x2)) + 
     &            muu3**2*(-8*p1k2**2*sqrtS**2 - 2*p1k2*sqrtS**4*(x1 - 2*(1 + x2)) + sqrtS**6*(-2 + x1 - 4*muu1**2*(-1 + x1 + x2))))*FL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(8*p1k2**2 + 2*p1k2*sqrtS**2*(1 + 8*muu3**2 - 2*x1 - 5*x2) + p1k1*(-8*p1k2 + 2*sqrtS**2*(3 + x2)) + 
     &            sqrtS**4*(-2 + x1 - 2*muu1**2*x1 + 2*muu1**2*x2 - 2*muu3**2*(4 + x1 + 3*x2) + x2*(x1 + 4*x2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(8*p1k2**2 + 2*p1k2*sqrtS**2*(1 + 8*muu3**2 - 2*x1 - 5*x2) + p1k1*(-8*p1k2 + 2*sqrtS**2*(3 + x2)) + 
     &            sqrtS**4*(-2 + x1 - 2*muu1**2*x1 + 2*muu1**2*x2 - 2*muu3**2*(4 + x1 + 3*x2) + x2*(x1 + 4*x2)))*FL(indj,i,tt,gen) + 
     &         (-4*muu3**4*sqrtS**6 + 2*p1k1*(8*p1k2**2 + 2*p1k2*sqrtS**2*(-1 + 2*muu3**2 - 3*x2) + sqrtS**4*x2*(1 + muu1**2 - muu3**2 + x2)) + 
     &            muu1**2*sqrtS**4*(-2*p1k2*(-2 + x1 + 2*x2) + sqrtS**2*x2*(-3 + x1 + 3*x2)) + 
     &            muu3**2*(-8*p1k2**2*sqrtS**2 - 2*p1k2*sqrtS**4*(x1 - 2*(1 + x2)) + sqrtS**6*(-2 + x1 - 4*muu1**2*(-1 + x1 + x2))))*FR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x2)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*(-2*p1k2 - sqrtS**2*(muu1**2 + muu3**2 - x2) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))  )

	endif	

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2d0*dreal( -(-2*CF*gs**2*NC*(dconjg(AL)*(2*dconjg(BR)*((muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + 
     &               p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + muu1**2*sqrtS**4*(1 + muu4**2 - x2))*ELc(indi,i,tt,gen) + 
     &            muu1*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*ER(indi,i,tt,gen))*
     &          (-2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + (2*p1k2 - sqrtS**2*x2)*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu3*ELc(indi,i,tt,gen)*(muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) + sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                FL(indj,i,tt,gen) + muu4*sqrtS**2*(-2*(sqrtS**2 + muu1**4*sqrtS**2 - 2*muu1**2*(p1k2 + sqrtS**2) + p1k2*x1) + 
     &                  muu2**2*(-4*p1k1 + sqrtS**2*(2 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) + 2*p1k1*x2 + 
     &                  sqrtS**2*((-1 + muu1**2)*(2*(muu3**2 + muu4**2) - x1) + 2*x2 - 4*muu1**2*x2))*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x2)*(-1 + x1 + x2) + 
     &                  muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*FL(indj,i,tt,gen) + 
     &               muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-1 - muu4**2*(1 + x1) + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*FR(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(2*dconjg(BL)*(muu1*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*muu3**2 - 2*muu4**2 + 2*x2))*ELc(indi,i,tt,gen) - 
     &            (muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + 
     &               muu1**2*sqrtS**4*(1 + muu4**2 - x2))*ER(indi,i,tt,gen))*((-2*p1k2 + sqrtS**2*x2)*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-1 - muu4**2*(1 + x1) + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x2)*(-1 + x1 + x2) + 
     &                  muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*FR(indj,i,tt,gen))
     &             + muu3*ER(indi,i,tt,gen)*(muu4*sqrtS**2*(-2*(sqrtS**2 + muu1**4*sqrtS**2 - 2*muu1**2*(p1k2 + sqrtS**2) + p1k2*x1) + 
     &                  muu2**2*(-4*p1k1 + sqrtS**2*(2 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) + 2*p1k1*x2 + 
     &                  sqrtS**2*((-1 + muu1**2)*(2*(muu3**2 + muu4**2) - x1) + 2*x2 - 4*muu1**2*x2))*FL(indj,i,tt,gen) + 
     &               muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + sqrtS**2*
     &                   (-2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) + sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmZ)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (-2*p1k2 - sqrtS**2*(muu2**2 + muu4**2 - x2) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

** TSchannel bremsstrahlung 32

        if (icha.eq.2) then
		
** TS-channel photon (for cha cha -> t tbar)
** relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(2*ApL*BpL*CF*gs**2*NC*(ER(indi,i,tt,gen)*((-8*p1k1**2*(2*p1k2 + sqrtS**2*(muu3**2 - x2)) + 
     &            2*p1k1*sqrtS**2*(2*p1k2*(-1 + 2*muu3**2 + x1) + sqrtS**2*(2*muu3**2*(-1 + x1) + x2 - (muu3**2 + x1)*x2 + muu1**2*(-2 + 2*x1 + x2))) + 
     &            sqrtS**4*(-4*muu3**4*sqrtS**2 + muu1**2*x1*(-2*p1k2 + sqrtS**2*(-1 + x1 + x2)) + muu3**2*(-2*p1k2*x1 + sqrtS**2*(-2 + 2*x1 + x2 - 4*muu1**2*(-1 + x1 + x2)))))*
     &          FL(indj,i,tt,gen) + muu1*muu3*sqrtS**2*(8*p1k1**2 - 2*p1k1*(4*p1k2 + sqrtS**2*(1 + 8*muu3**2 - x1 - 4*x2)) + 
     &            sqrtS**2*(2*p1k2*(-3 + x1) + sqrtS**2*(-2 + x1 + 2*muu1**2*x1 + x1**2 + 2*muu3**2*(-4 + x1 - x2) - 2*(-2 + muu1**2 + x1)*x2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(8*p1k1**2 - 2*p1k1*(4*p1k2 + sqrtS**2*(1 + 8*muu3**2 - x1 - 4*x2)) + 
     &            sqrtS**2*(2*p1k2*(-3 + x1) + sqrtS**2*(-2 + x1 + 2*muu1**2*x1 + x1**2 + 2*muu3**2*(-4 + x1 - x2) - 2*(-2 + muu1**2 + x1)*x2)))*FL(indj,i,tt,gen) + 
     &         (-8*p1k1**2*(2*p1k2 + sqrtS**2*(muu3**2 - x2)) + 2*p1k1*sqrtS**2*
     &             (2*p1k2*(-1 + 2*muu3**2 + x1) + sqrtS**2*(2*muu3**2*(-1 + x1) + x2 - (muu3**2 + x1)*x2 + muu1**2*(-2 + 2*x1 + x2))) + 
     &            sqrtS**4*(-4*muu3**4*sqrtS**2 + muu1**2*x1*(-2*p1k2 + sqrtS**2*(-1 + x1 + x2)) + muu3**2*(-2*p1k2*x1 + sqrtS**2*(-2 + 2*x1 + x2 - 4*muu1**2*(-1 + x1 + x2)))))*
     &          FR(indj,i,tt,gen))))/(sqrtS**4*(-1 + x1)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu1**2 + muu3**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))    )

	endif	

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  - (-2*CF*gs**2*NC*(dconjg(AL)*(2*dconjg(BR)*(p1k1*ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen))*
     &          (2*muu2*muu4*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*FL(indj,i,tt,gen) - 
     &            (4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*(-2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + 
     &                  sqrtS**2*(-2*muu4**2*(-1 + muu1**2 + x1) + 2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + muu4**2 + x1)*x2)))*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu3*ELc(indi,i,tt,gen)*(muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*FL(indj,i,tt,gen) + 
     &               muu4*sqrtS**2*(-2*muu2**4*sqrtS**2 + 2*muu1**2*(2*p1k2 + sqrtS**2) - 2*(sqrtS**2 + p1k2*x1) + 
     &                  muu2**2*(-4*p1k1 + sqrtS**2*(2*(2 + muu1**2 + muu3**2 + muu4**2 - x1) - x2)) + 2*p1k1*x2 + sqrtS**2*(2*(-1 + muu1**2)*(muu3**2 + muu4**2 - x1) + x2 - 3*muu1**2*x2))*
     &                FR(indj,i,tt,gen)) + muu1*ER(indi,i,tt,gen)*(muu2*sqrtS**2*
     &                (muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + 
     &                  muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x1)*(-1 + x1 + x2))*FL(indj,i,tt,gen) + 
     &               muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-((1 + muu4**2 - x1)*(1 + x1)) + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*FR(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(2*dconjg(BL)*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen))*
     &          (-((4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*(-2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + 
     &                    sqrtS**2*(-2*muu4**2*(-1 + muu1**2 + x1) + 2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + muu4**2 + x1)*x2)))*FL(indj,i,tt,gen)) + 
     &            2*muu2*muu4*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-((1 + muu4**2 - x1)*(1 + x1)) + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + 
     &                  muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x1)*(-1 + x1 + x2))*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*(muu4*sqrtS**2*(-2*muu2**4*sqrtS**2 + 2*muu1**2*(2*p1k2 + sqrtS**2) - 2*(sqrtS**2 + p1k2*x1) + 
     &                  muu2**2*(-4*p1k1 + sqrtS**2*(2*(2 + muu1**2 + muu3**2 + muu4**2 - x1) - x2)) + 2*p1k1*x2 + sqrtS**2*(2*(-1 + muu1**2)*(muu3**2 + muu4**2 - x1) + x2 - 3*muu1**2*x2))*
     &                FL(indj,i,tt,gen) + muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmZ)*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

        if (cutoff.lt.zerocheck) then

** TSchannel dipole with Z 312

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal(  (-4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (dconjg(AR)*(-4*dconjg(BL)*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + 
     &            (((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &               (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/(4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &             ER(indi,i,tt,gen))*(-((-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                  sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                     x2*((-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**4*x2 + 
     &                        (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FL(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2)) + muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(-(muu4*(-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                     sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                        x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                           (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FL(indj,i,tt,gen))/
     &                (2.*sqrtS**4*(4*muu4**2 - x2**2)) + muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FL(indj,i,tt,gen) - 
     &               2*muu2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &                  (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &                  ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/(4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &                FR(indj,i,tt,gen)))) + dconjg(AL)*(-4*dconjg(BR)*((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/(4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &             ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen))*
     &          (muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) - ((-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                    x2*((-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**4*x2 + 
     &                       (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FR(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2))) + sqrtS**2*dconjg(BL)*
     &          (muu3*ELc(indi,i,tt,gen)*(-2*muu2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &                  (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &                  ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/(4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*
     &                FL(indj,i,tt,gen) + (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*FL(indj,i,tt,gen) - 
     &               (muu4*(-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                    sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                       x2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**4*x2 + 
     &                          (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*FR(indj,i,tt,gen))/
     &                (2.*sqrtS**4*(4*muu4**2 - x2**2)))))))/
     &  (sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*dconjg(DenSmZ)*((muu1**2 + muu3**2)*sqrtS**2 - 
     &      2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &         (p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x2)/(4.*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))) - 
     &      (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

** TSchannel dipole with Z 321

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal( (-4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (3*muu3**4 - muu4**4 + muu3**2*(4 + 2*muu4**2 - 5*x1 - x2) + muu4**2*(-x1 + x2) + (-1 + x1)*(-1 + 2*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (dconjg(AR)*(-4*dconjg(BL)*(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + 
     &            (((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &               (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &             ER(indi,i,tt,gen))*(((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                 sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                    x1*((-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**4*x1 + 
     &                       (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FL(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu3**2 - x1**2)) + muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu1*ELc(indi,i,tt,gen)*(-(muu4*(4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                     sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                        x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &                           (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FL(indj,i,tt,gen))/
     &                (2.*sqrtS**4*(4*muu3**2 - x1**2)) + muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FL(indj,i,tt,gen) - 
     &               2*muu2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &                  (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &                FR(indj,i,tt,gen)))) + dconjg(AL)*(-4*dconjg(BR)*((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - 
     &               ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &               ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &             ELc(indi,i,tt,gen) + muu1*muu3*sqrtS**2*ER(indi,i,tt,gen))*
     &          (muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + ((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                 sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                    x1*((-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**4*x1 + 
     &                       (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FR(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu3**2 - x1**2))) + sqrtS**2*dconjg(BL)*
     &          (muu3*ELc(indi,i,tt,gen)*(-2*muu2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. - 
     &                  (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &                FL(indj,i,tt,gen) + (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*FL(indj,i,tt,gen) - 
     &               (muu4*(4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                    sqrtS**2*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**4 - (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                       x1*((-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**4*x1 + 
     &                          (1 + muu1**2 - muu2**2)*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*FR(indj,i,tt,gen))/
     &                (2.*sqrtS**4*(4*muu3**2 - x1**2)))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmZ)*((muu1**2 + muu3**2)*sqrtS**2 - 
     &      2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))) - 
     &      (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))) )

        if(icha.eq.2) then

** TSchannel dipole with G 312

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal( (4*ApL*BpL*CF*GS**2*NC*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))*(-2*muu3**2*
     &       (4*Sqrt(1 - 4*muu3**2) - 8*Sqrt(1 - 4*muu3**2)*x2 + 3*Sqrt(1 - 4*muu3**2)*x2**2 + Sqrt(1 - 4*muu3**2)*x1*(-2 + 3*x2) - 4*Sqrt(-4*muu3**2 + x2**2) + 
     &         4*x2*Sqrt(-4*muu3**2 + x2**2)) + (-1 + x2)*(Sqrt(1 - 4*muu3**2)*x1**2 + 3*Sqrt(1 - 4*muu3**2)*x1*(-1 + x2) - 5*Sqrt(1 - 4*muu3**2)*x2 + 2*Sqrt(1 - 4*muu3**2)*x2**2 + 
     &         2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x2**2))) + 4*muu3**4*
     &       (Sqrt(1 - 4*muu3**2)*x1 - 2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x2**2)) + x2*(Sqrt(1 - 4*muu3**2) + 2*Sqrt(-4*muu3**2 + x2**2))))*
     &    (ER(indi,i,tt,gen)*((8*(1 - 4*muu3**2)*p1k2**2*sqrtS**2 - 2*(1 + 4*muu1**2)*muu3**2*sqrtS**6 - 8*muu3**4*sqrtS**6 + 
     &            sqrtS**2*x2*(sqrtS**4*(x2 + 2*muu1**2*x2) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))) - 
     &            4*p1k2*(sqrtS**4*(x2 - 4*muu3**2*x2) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))))*FL(indj,i,tt,gen) + 
     &         2*muu1*muu3*(-12*muu3**2*sqrtS**6 + 3*sqrtS**6*x2**2 - 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*((-24*muu1*muu3**3*sqrtS**6 + 6*muu1*muu3*sqrtS**6*x2**2 - 8*muu1*muu3*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            2*muu1*muu3*sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)))*FL(indj,i,tt,gen) + 
     &         (8*(1 - 4*muu3**2)*p1k2**2*sqrtS**2 - 2*(1 + 4*muu1**2)*muu3**2*sqrtS**6 - 8*muu3**4*sqrtS**6 + 
     &            sqrtS**2*x2*(sqrtS**4*(x2 + 2*muu1**2*x2) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))) - 
     &            4*p1k2*(sqrtS**4*(x2 - 4*muu3**2*x2) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))))*FR(indj,i,tt,gen))))/
     &  ((-1 + 2*muu3**2)*sqrtS**6*(-1 + x2)**2*(-2 + x1 + x2)*(-4*muu3**2 + x2**2)**1.5*
     &    (-4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4) + sqrtS**2*(Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2 + (1 - 2*muu1**2 - 2*muu3**2)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))) + 
     &      2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2))*(MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))) )

** TSchannel dipole with G 321

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal((4*ApL*BpL*CF*GS**2*NC*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*(4*muu3**4*
     &       (-2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x1**2)) + x1*(Sqrt(1 - 4*muu3**2) + 2*Sqrt(-4*muu3**2 + x1**2)) + Sqrt(1 - 4*muu3**2)*x2) + 
     &      (-1 + x1)*(2*Sqrt(1 - 4*muu3**2)*x1**2 + 2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x1**2)) - 3*Sqrt(1 - 4*muu3**2)*x2 + Sqrt(1 - 4*muu3**2)*x2**2 + 
     &         Sqrt(1 - 4*muu3**2)*x1*(-5 + 3*x2)) - 2*muu3**2*(4*Sqrt(1 - 4*muu3**2) + 3*Sqrt(1 - 4*muu3**2)*x1**2 - 4*Sqrt(-4*muu3**2 + x1**2) - 2*Sqrt(1 - 4*muu3**2)*x2 + 
     &         x1*(-8*Sqrt(1 - 4*muu3**2) + 4*Sqrt(-4*muu3**2 + x1**2) + 3*Sqrt(1 - 4*muu3**2)*x2)))*
     &    (ER(indi,i,tt,gen)*((8*(-1 + 4*muu3**2)*p1k1**2*sqrtS**2 + 2*(1 + 4*muu1**2)*muu3**2*sqrtS**6 + 8*muu3**4*sqrtS**6 - 
     &            4*p1k1*((-1 + 4*muu3**2)*sqrtS**4*x1 + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &            sqrtS**2*x1*(-(sqrtS**4*(x1 + 2*muu1**2*x1)) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*FL(indj,i,tt,gen) + 
     &         2*muu1*muu3*(12*muu3**2*sqrtS**6 - 3*sqrtS**6*x1**2 - 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))*FR(indj,i,tt,gen)) + 
     &      ELc(indi,i,tt,gen)*(2*muu1*muu3*(12*muu3**2*sqrtS**6 - 3*sqrtS**6*x1**2 - 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))*FL(indj,i,tt,gen) + 
     &         (8*(-1 + 4*muu3**2)*p1k1**2*sqrtS**2 + 2*(1 + 4*muu1**2)*muu3**2*sqrtS**6 + 8*muu3**4*sqrtS**6 - 
     &            4*p1k1*((-1 + 4*muu3**2)*sqrtS**4*x1 + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &            sqrtS**2*x1*(-(sqrtS**4*(x1 + 2*muu1**2*x1)) + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*FR(indj,i,tt,gen))))/
     &  ((-1 + 2*muu3**2)*sqrtS**6*(-1 + x1)**2*(-4*muu3**2 + x1**2)**1.5*(-2 + x1 + x2)*
     &    (-4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4) + sqrtS**2*(Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1 + (-1 + 2*muu1**2 + 2*muu3**2)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) - 
     &      2*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*(MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))))

         endif
 
        endif


	  do j = jmin, jmax

** TSchannel bremsstrahlung 13

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(   -(CF*gs**2*NC*(dconjg(CL(j))*(2*dconjg(DR(j))*(muu3*(1 + muu3**2 + muu4**2 - x2)*ELc(indi,i,tt,gen)*
     &             (-((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FL(indj,i,tt,gen)) + 2*muu2*p1k2*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(-(muu4*((1 - muu1**2 + muu2**2)*sqrtS**2 + 2*p1k2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*FL(indj,i,tt,gen)) + 
     &               muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*FR(indj,i,tt,gen))) + 
     &         dconjg(DL(j))*(-2*muu1*muu3*(1 + muu3**2 + muu4**2 - x2)*ER(indi,i,tt,gen)*((-2*p1k2 + sqrtS**2*x2)*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*((2*p1k2*(-(muu3**2*(-6 + x1)) + (-2 + x1)*(1 + muu4**2 - x2)) + 2*p1k1*x2*(-1 + muu3**2 - muu4**2 + x2) + 
     &                  sqrtS**2*((-1 + muu2)*(1 + muu2)*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 6*x2) + (1 + muu4**2 - x2)*(-1 + muu4**2 + x1 - x2)) + 
     &                     muu1**2*(muu3**4 - muu3**2*(6*muu4**2 + x1) + (1 + muu4**2 - x2)*(-1 + muu4**2 + x1 + x2))))*FL(indj,i,tt,gen) - 
     &               2*muu2*muu4*(-2*p1k1*(-1 + muu3**2 - muu4**2 + x2) - (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*FR(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(-2*dconjg(DL(j))*(muu3*(1 + muu3**2 + muu4**2 - x2)*ER(indi,i,tt,gen)*
     &             (-2*muu2*p1k2*FL(indj,i,tt,gen) + (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*FR(indj,i,tt,gen)) + 
     &            muu1*ELc(indi,i,tt,gen)*(-(muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*FL(indj,i,tt,gen)) + 
     &               muu4*((1 - muu1**2 + muu2**2)*sqrtS**2 + 2*p1k2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                  sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*FR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(-2*muu1*muu3*(1 + muu3**2 + muu4**2 - x2)*ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*FL(indj,i,tt,gen) + (-2*p1k2 + sqrtS**2*x2)*FR(indj,i,tt,gen)) + 
     &            ER(indi,i,tt,gen)*(-2*muu2*muu4*(-2*p1k1*(-1 + muu3**2 - muu4**2 + x2) - (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*FL(indj,i,tt,gen) + 
     &               (2*p1k2*(-(muu3**2*(-6 + x1)) + (-2 + x1)*(1 + muu4**2 - x2)) + 2*p1k1*x2*(-1 + muu3**2 - muu4**2 + x2) + 
     &                  sqrtS**2*((-1 + muu2)*(1 + muu2)*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 6*x2) + (1 + muu4**2 - x2)*(-1 + muu4**2 + x1 - x2)) + 
     &                     muu1**2*(muu3**4 - muu3**2*(6*muu4**2 + x1) + (1 + muu4**2 - x2)*(-1 + muu4**2 + x1 + x2))))*FR(indj,i,tt,gen))))))/
     &  ((-1 + muu3**2 - muu4**2 + x2)**2*dconjg(DenSmH(j))*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))    )

** TSchannel bremsstrahlung 14

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(   ((CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(muu1*muu3*ER(indi,i,tt,gen)*
     &               ((-2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*x2)) + muu4**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x2)) + 2*p1k1*x2 - 
     &                    (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FR(indj,i,tt,gen)) + 
     &              ELc(indi,i,tt,gen)*(((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 - 4*p1k2*x1 - 2*(-1 + muu2**2)*sqrtS**2*x1*x2 + 
     &                    muu4**2*(2*p1k2*x1 - 2*p1k1*(-2 + x2) + sqrtS**2*(muu1**2*(1 - 6*muu3**2 + x1 - x2) - (-1 + muu2**2)*(-1 + 6*muu3**2 - 3*x1 + x2))) - 
     &                    muu3**2*(2*p1k2*(-6 + x1) - 2*p1k1*x2 + sqrtS**2*((-1 + muu2**2)*(-1 + x1 - 5*x2) + muu1**2*(-1 + x1 + x2))))*FL(indj,i,tt,gen) + 
     &                 muu2*muu4*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 + muu3**2*(4*p1k1 + 12*p1k2 - 5*(1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                    muu4**2*(-4*(p1k1 + p1k2) + (1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 2*p1k1*x2 + 
     &                    (4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*x2)*FR(indj,i,tt,gen))) + 
     &           dconjg(DR(j))*(-(muu1*ER(indi,i,tt,gen)*(muu4*(4*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 + (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2)*x1 + 
     &                      muu4**2*(4*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) - (2*p1k1 + 4*p1k2 + (3 - muu1**2 + muu2**2)*sqrtS**2)*x2 + 2*sqrtS**2*x2**2 + 
     &                      muu3**2*(-4*p1k1 - 12*p1k2 + sqrtS**2*(-5 + 5*muu1**2 - 5*muu2**2 + 2*x1 + 6*x2)))*FL(indj,i,tt,gen) - 
     &                   2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*FR(indj,i,tt,gen))) + 
     &              muu3*ELc(indi,i,tt,gen)*(-2*muu4*(-2*p1k1 - 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 + x1 + x2 - muu2**2*(x1 + x2)))*
     &                  FL(indj,i,tt,gen) + muu2*(muu3**2*(4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu4**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 2*p1k1*x2 - 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FR(indj,i,tt,gen)))) + 
     &        dconjg(CR(j))*(dconjg(DL(j))*(muu3*ER(indi,i,tt,gen)*(muu2*(muu3**2*(4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu4**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 
     &                    2*p1k2*x1 + 2*p1k1*x2 - (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FL(indj,i,tt,gen) - 
     &                 2*muu4*(-2*p1k1 - 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 + x1 + x2 - muu2**2*(x1 + x2)))*FR(indj,i,tt,gen)) + 
     &              muu1*ELc(indi,i,tt,gen)*(2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*FL(indj,i,tt,gen) - 
     &                 muu4*(4*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 + (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2)*x1 + 
     &                    muu4**2*(4*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) - (2*p1k1 + 4*p1k2 + (3 - muu1**2 + muu2**2)*sqrtS**2)*x2 + 2*sqrtS**2*x2**2 + 
     &                    muu3**2*(-4*p1k1 - 12*p1k2 + sqrtS**2*(-5 + 5*muu1**2 - 5*muu2**2 + 2*x1 + 6*x2)))*FR(indj,i,tt,gen))) + 
     &           dconjg(DR(j))*(muu1*muu3*ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FL(indj,i,tt,gen) + 
     &                 (-2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*x2)) + muu4**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x2)) + 2*p1k1*x2 - 
     &                    (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FR(indj,i,tt,gen)) + 
     &              ER(indi,i,tt,gen)*(muu2*muu4*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 + muu3**2*(4*p1k1 + 12*p1k2 - 5*(1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                    muu4**2*(-4*(p1k1 + p1k2) + (1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 2*p1k1*x2 + 
     &                    (4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*x2)*FL(indj,i,tt,gen) + 
     &                 ((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 - 4*p1k2*x1 - 2*(-1 + muu2**2)*sqrtS**2*x1*x2 + 
     &                    muu4**2*(2*p1k2*x1 - 2*p1k1*(-2 + x2) + sqrtS**2*(muu1**2*(1 - 6*muu3**2 + x1 - x2) - (-1 + muu2**2)*(-1 + 6*muu3**2 - 3*x1 + x2))) - 
     &                    muu3**2*(2*p1k2*(-6 + x1) - 2*p1k1*x2 + sqrtS**2*((-1 + muu2**2)*(-1 + x1 - 5*x2) + muu1**2*(-1 + x1 + x2))))*FR(indj,i,tt,gen))))))/
     &    ((1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmH(j))*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))))   )

** TSchannel bremsstrahlung 23

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal( ((CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DR(j))*(muu1*ER(indi,i,tt,gen)*
     &               (muu4*(-(muu3**2*(4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1))) - 2*p1k2*x1 + muu4**2*(-4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*x1)) + 2*p1k1*x2 + 
     &                    (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FL(indj,i,tt,gen) + 
     &                 2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*FR(indj,i,tt,gen)) + 
     &              muu3*ELc(indi,i,tt,gen)*(-2*muu4*(2*p1k1 + 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 - muu1**2*(x1 + x2)))*
     &                  FL(indj,i,tt,gen) + muu2*(-4*muu4**2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2 + 5*muu4**2*sqrtS**2 + 5*muu1**2*muu4**2*sqrtS**2 - 5*muu2**2*muu4**2*sqrtS**2 + 
     &                    muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + (-1 - muu1**2 + muu2**2)*sqrtS**2*x2 + 
     &                    2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2))*FR(indj,i,tt,gen))) + 
     &           dconjg(DL(j))*(muu1*muu3*ER(indi,i,tt,gen)*((-4*muu4**2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 - 5*muu4**2*sqrtS**2 + 5*muu1**2*muu4**2*sqrtS**2 - 
     &                    5*muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 3*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 6*muu4**2*sqrtS**2*x1 + 2*sqrtS**2*x1**2 + 
     &                    muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) + (1 - muu1**2 + muu2**2 + 2*muu4**2)*sqrtS**2*x2 + 2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2)
     &                    )*FL(indj,i,tt,gen) + 2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FR(indj,i,tt,gen)) + 
     &              ELc(indi,i,tt,gen)*(((-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + muu3**2*((-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 - x1) + 2*p1k2*(-2 + x1)) + 
     &                    (-2*(-2 + muu3**2)*p1k1 + sqrtS**2*((-1 + 3*muu1**2 + muu2**2)*muu3**2 - 2*muu1**2*x1))*x2 + 
     &                    muu4**2*(-2*p1k2*x1 + 2*p1k1*(-6 + x2) + sqrtS**2*(-(muu1**2*(-1 + 6*muu3**2 - 5*x1 + x2)) - (-1 + muu2**2)*(-1 + 6*muu3**2 + x1 + x2))))*FL(indj,i,tt,gen) + 
     &                 muu2*muu4*(muu3**2*(-4*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu4**2*(-4*p1k1 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 2*p1k1*x2 + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FR(indj,i,tt,gen)))) + 
     &        dconjg(CR(j))*(dconjg(DL(j))*(muu1*ELc(indi,i,tt,gen)*(2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*
     &                  FL(indj,i,tt,gen) + muu4*(-(muu3**2*(4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1))) - 2*p1k2*x1 + muu4**2*(-4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*x1)) + 
     &                    2*p1k1*x2 + (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FR(indj,i,tt,gen)) + 
     &              muu3*ER(indi,i,tt,gen)*(muu2*(-4*muu4**2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2 + 5*muu4**2*sqrtS**2 + 5*muu1**2*muu4**2*sqrtS**2 - 5*muu2**2*muu4**2*sqrtS**2 + 
     &                    muu3**2*(4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + (-1 - muu1**2 + muu2**2)*sqrtS**2*x2 + 
     &                    2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2))*FL(indj,i,tt,gen) - 
     &                 2*muu4*(2*p1k1 + 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 - muu1**2*(x1 + x2)))*FR(indj,i,tt,gen))) + 
     &           dconjg(DR(j))*(muu1*muu3*ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*FL(indj,i,tt,gen) + 4*p1k1*FR(indj,i,tt,gen) + 
     &                 (-2*(2*p1k1 + p1k2)*x1 + muu3**2*(4*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) + 2*p1k1*x2 + 
     &                    sqrtS**2*((-1 + x1)*(-1 + muu1**2 - muu2**2 + 2*x1) + (1 - muu1**2 + muu2**2)*x2) + 
     &                    muu4**2*(-12*p1k1 - 4*p1k2 + sqrtS**2*(-5 + 5*muu1**2 - 5*muu2**2 + 6*x1 + 2*x2)))*FR(indj,i,tt,gen)) + 
     &              ER(indi,i,tt,gen)*(muu2*muu4*(muu3**2*(-4*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu4**2*(-4*p1k1 + (-1 - muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 2*p1k1*x2 + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*FL(indj,i,tt,gen) + 
     &                 ((-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + muu3**2*((-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 - x1) + 2*p1k2*(-2 + x1)) + 
     &                    (-2*(-2 + muu3**2)*p1k1 + sqrtS**2*((-1 + 3*muu1**2 + muu2**2)*muu3**2 - 2*muu1**2*x1))*x2 + 
     &                    muu4**2*(-2*p1k2*x1 + 2*p1k1*(-6 + x2) + sqrtS**2*(-(muu1**2*(-1 + 6*muu3**2 - 5*x1 + x2)) - (-1 + muu2**2)*(-1 + 6*muu3**2 + x1 + x2))))*FR(indj,i,tt,gen))))))/
     &    ((1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmH(j))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))))   )

** TSchannel bremsstrahlung 24

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(  -(CF*gs**2*NC*(dconjg(CR(j))*(dconjg(DR(j))*(-4*muu2*muu4*(1 + muu3**2 + muu4**2 - x1)
     &            *(muu1*muu3*sqrtS**2*ELc(indi,i,tt,gen) + p1k1*ER(indi,i,tt,gen))*FL(indj,i,tt,gen) + 
     &            (-2*muu1*muu3*(2*p1k2*(-1 - muu3**2 + muu4**2 + x1) + 2*p1k1*(-1 - muu3**2 + 3*muu4**2 + x1) + 
     &                  sqrtS**2*((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1) + (1 + muu3**2 - muu4**2 - x1)*x2))*ELc(indi,i,tt,gen) + 
     &               (2*p1k2*(1 + muu3**2 - muu4**2 - x1)*x1 - 2*p1k1*(-2 + 6*muu4**2 + 2*x1 + muu3**2*(-2 + x2) + x2 - (muu4**2 + x1)*x2) + 
     &                  sqrtS**2*((-1 + muu2**2)*(muu3**4 + (-1 + muu4**2 + x1)*(1 + muu4**2 - x1 - x2) + muu3**2*(-6*muu4**2 + x2)) + 
     &                     muu1**2*(-1 + muu3**4 + muu4**4 + 6*muu4**2*x1 + x1**2 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - 2*x1 + x2))))*ER(indi,i,tt,gen))*FR(indj,i,tt,gen)) - 
     &         2*dconjg(DL(j))*(muu3*ER(indi,i,tt,gen)*(muu2*((1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(-1 - muu3**2 + muu4**2 + x1) + 
     &                  2*p1k1*(-1 - muu3**2 + 3*muu4**2 + x1))*FL(indj,i,tt,gen) + (-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FR(indj,i,tt,gen)) + 
     &            muu1*ELc(indi,i,tt,gen)*(-(muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*FL(indj,i,tt,gen)) + 
     &               muu4*(1 + muu3**2 + muu4**2 - x1)*(2*p1k1 - sqrtS**2*x1)*FR(indj,i,tt,gen)))) + 
     &      dconjg(CL(j))*(dconjg(DL(j))*(ELc(indi,i,tt,gen)*((2*p1k2*(1 + muu3**2 - muu4**2 - x1)*x1 - 2*p1k1*(-2 + 6*muu4**2 + 2*x1 + muu3**2*(-2 + x2) + x2 - (muu4**2 + x1)*x2) + 
     &                  sqrtS**2*((-1 + muu2**2)*(muu3**4 + (-1 + muu4**2 + x1)*(1 + muu4**2 - x1 - x2) + muu3**2*(-6*muu4**2 + x2)) + 
     &                     muu1**2*(-1 + muu3**4 + muu4**4 + 6*muu4**2*x1 + x1**2 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - 2*x1 + x2))))*FL(indj,i,tt,gen) - 
     &               4*muu2*muu4*p1k1*(1 + muu3**2 + muu4**2 - x1)*FR(indj,i,tt,gen)) + 
     &            2*muu1*muu3*ER(indi,i,tt,gen)*((2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                  sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*FL(indj,i,tt,gen) - 
     &               2*muu2*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FR(indj,i,tt,gen))) + 
     &         2*dconjg(DR(j))*(muu3*ELc(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*FL(indj,i,tt,gen)) + 
     &               muu2*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + 3*muu4**2 + x1))*FR(indj,i,tt,gen))
     &             + muu1*ER(indi,i,tt,gen)*(muu4*(1 + muu3**2 + muu4**2 - x1)*(-2*p1k1 + sqrtS**2*x1)*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*FR(indj,i,tt,gen))))))/
     &  ((-1 - muu3**2 + muu4**2 + x1)**2*dconjg(DenSmH(j))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))  )

** TSchannel bremsstrahlung 33

	ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS + 2*dreal(   -(CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(2*muu1*muu3*ER(indi,i,tt,gen)*
     &             ((4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*(-2*p1k2*(-1 + muu3**2 + muu4**2 + 2*x1 - x2) + 
     &                     sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1) + (-1 + muu3**2 + 2*x1 - x2)*x2 + muu4**2*(-2*(-1 + muu1**2 + x1) + x2))))*FL(indj,i,tt,gen) + 
     &               2*muu2*muu4*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*((2*p1k1*(8*p1k2 + sqrtS**2*((-3 - 2*muu1**2 + 4*muu2**2 + muu3**2 - x2)*x2 + muu4**2*(4 - 4*muu2**2 + x2))) + 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 
     &                     muu3**2*(p1k2*(8*muu1**2 - 2*(4 + x1)) + sqrtS**2*(-2*muu1**4 + (-1 + muu2**2)*(2 + 6*muu4**2 - x1 - 6*x2) + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - x1 - 6*x2))) - 
     &                     2*p1k2*x1*(1 + 2*muu1**2 + muu4**2 - x2) + sqrtS**2*
     &                      ((-1 + muu2**2)*(1 + muu4**2 - x2)*(1 + muu4**2 - x1 - x2) + 2*muu1**4*(-1 - muu4**2 + x1 + x2) + 
     &                        muu1**2*(3 + muu4**4 - 3*x1 + muu4**2*(4 - 5*x1 - 2*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 4*x2 + 5*x1*x2 + x2**2))))*FL(indj,i,tt,gen) + 
     &               4*muu2*muu4*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + 
     &                  muu1**2*sqrtS**4*(1 + muu4**2 - x2))*FR(indj,i,tt,gen))) + 
     &         2*dconjg(DR(j))*(muu1*ER(indi,i,tt,gen)*(muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-1 - muu4**2*(1 + x1) + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*FL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x2)*(-1 + x1 + x2) + 
     &                  muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*FR(indj,i,tt,gen))
     &             + muu3*ELc(indi,i,tt,gen)*(muu4*sqrtS**2*(-2*(sqrtS**2 + muu1**4*sqrtS**2 - 2*muu1**2*(p1k2 + sqrtS**2) + p1k2*x1) + 
     &                  muu2**2*(-4*p1k1 + sqrtS**2*(2 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) + 2*p1k1*x2 + 
     &                  sqrtS**2*((-1 + muu1**2)*(2*(muu3**2 + muu4**2) - x1) + 2*x2 - 4*muu1**2*x2))*FL(indj,i,tt,gen) + 
     &               muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + sqrtS**2*
     &                   (-2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) + sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                FR(indj,i,tt,gen)))) + dconjg(CR(j))*(-2*dconjg(DL(j))*
     &          (muu3*ER(indi,i,tt,gen)*(muu2*(p1k1*(-8*p1k2 + 2*sqrtS**2*(-2*muu4**2 + x2)) + 
     &                  sqrtS**2*(2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) - sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                FL(indj,i,tt,gen) + muu4*sqrtS**2*(-((-1 + muu1**2)*sqrtS**2*(2*(muu3**2 + muu4**2) - x1)) + 2*(sqrtS**2 + muu1**4*sqrtS**2 - 2*muu1**2*(p1k2 + sqrtS**2) + p1k2*x1) - 
     &                  2*(p1k1 + (1 - 2*muu1**2)*sqrtS**2)*x2 + muu2**2*(4*p1k1 + sqrtS**2*(-2 - 2*muu1**2 - 2*muu3**2 - 2*muu4**2 + x1 + 2*x2)))*FR(indj,i,tt,gen)) + 
     &            muu1*ELc(indi,i,tt,gen)*(-(muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x2)*(-1 + x1 + x2) + 
     &                    muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*FL(indj,i,tt,gen))
     &                 + muu4*(-2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) - 
     &                     sqrtS**2*(-1 - muu4**2*(1 + x1) + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*FR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(2*muu1*muu3*ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*FL(indj,i,tt,gen) + 
     &               (4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*(-2*p1k2*(-1 + muu3**2 + muu4**2 + 2*x1 - x2) + 
     &                     sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1) + (-1 + muu3**2 + 2*x1 - x2)*x2 + muu4**2*(-2*(-1 + muu1**2 + x1) + x2))))*FR(indj,i,tt,gen)) + 
     &            ER(indi,i,tt,gen)*(4*muu2*muu4*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + 
     &                  muu1**2*sqrtS**4*(1 + muu4**2 - x2))*FL(indj,i,tt,gen) + 
     &               (2*p1k1*(8*p1k2 + sqrtS**2*((-3 - 2*muu1**2 + 4*muu2**2 + muu3**2 - x2)*x2 + muu4**2*(4 - 4*muu2**2 + x2))) + 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 
     &                     muu3**2*(p1k2*(8*muu1**2 - 2*(4 + x1)) + sqrtS**2*(-2*muu1**4 + (-1 + muu2**2)*(2 + 6*muu4**2 - x1 - 6*x2) + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - x1 - 6*x2))) - 
     &                     2*p1k2*x1*(1 + 2*muu1**2 + muu4**2 - x2) + sqrtS**2*
     &                      ((-1 + muu2**2)*(1 + muu4**2 - x2)*(1 + muu4**2 - x1 - x2) + 2*muu1**4*(-1 - muu4**2 + x1 + x2) + 
     &                        muu1**2*(3 + muu4**4 - 3*x1 + muu4**2*(4 - 5*x1 - 2*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 4*x2 + 5*x1*x2 + x2**2))))*FR(indj,i,tt,gen))))))/
     &  (2.*(-1 + muu3**2 - muu4**2 + x2)*dconjg(DenSmH(j))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))   )

** TSchannel bremsstrahlung 34

	ChiChi2QQ_RadiationTS =  ChiChi2QQ_RadiationTS + 2*dreal(  (CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(ELc(indi,i,tt,gen)
     &                   *((2*p1k1*(8*p1k2 + sqrtS**2*((-3 + 2*muu2**2 + muu3**2 - x1)*x2 + muu4**2*(4 - 4*muu2**2 + x2))) + 
     &                   sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-1 - 4*muu1**2 + 2*muu2**2 - muu4**2 + x1) + 
     &                      sqrtS**2*(-((-1 + muu2**2)*(-1 + 2*muu2**2 - muu4**2 + x1)*(1 + muu4**2 - x1 - x2)) + 
     &                         muu1**2*(muu4**4 + (-1 + x1)**2 + 2*muu2**2*(1 + muu4**2 - x1 - x2) + (-1 + 5*x1)*x2 - muu4**2*(-2 + 6*x1 + x2))) + 
     &                      muu3**2*(p1k2*(8*muu1**2 - 2*(4 + x1)) + sqrtS**2*(muu1**2*(2 + 2*muu2**2 + 6*muu4**2 - 2*x1 - 5*x2) - (-1 + muu2**2)*(2*(-1 + muu2**2 - 3*muu4**2 + x1) + 5*x2)))
     &                      ))*FL(indj,i,tt,gen) + 4*muu2*muu4*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 + muu3**2 + muu4**2 - x1)) + 
     &                   muu1**2*sqrtS**4*(1 + muu4**2 - x2))*FR(indj,i,tt,gen)) + 
     &             2*muu1*muu3*ER(indi,i,tt,gen)*((4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + 
     &                   sqrtS**2*(-2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + 
     &                      sqrtS**2*(-2*muu4**2*(-1 + muu1**2 + x1) + 2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + muu4**2 + x1)*x2)))*FL(indj,i,tt,gen) + 
     &                2*muu2*muu4*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*FR(indj,i,tt,gen))) + 
     &          2*dconjg(DR(j))*(muu1*ER(indi,i,tt,gen)*(muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                   sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                      sqrtS**2*(-((1 + muu4**2 - x1)*(1 + x1)) + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*FL(indj,i,tt,gen) + 
     &                muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + 
     &                   muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x1)*(-1 + x1 + x2))*FR(indj,i,tt,gen)) + 
     &             muu3*ELc(indi,i,tt,gen)*(muu4*sqrtS**2*(-2*muu2**4*sqrtS**2 + 2*muu1**2*(2*p1k2 + sqrtS**2) - 2*(sqrtS**2 + p1k2*x1) + 
     &                   muu2**2*(-4*p1k1 + sqrtS**2*(2*(2 + muu1**2 + muu3**2 + muu4**2 - x1) - x2)) + 2*p1k1*x2 + sqrtS**2*(2*(-1 + muu1**2)*(muu3**2 + muu4**2 - x1) + x2 - 3*muu1**2*x2))*
     &                 FL(indj,i,tt,gen) + muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                   sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                      sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*FR(indj,i,tt,gen)))) + 
     &       dconjg(CR(j))*(2*dconjg(DL(j))*(muu3*ER(indi,i,tt,gen)*(muu2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                   sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                      sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*FL(indj,i,tt,gen) + 
     &                muu4*sqrtS**2*(-2*muu2**4*sqrtS**2 + 2*muu1**2*(2*p1k2 + sqrtS**2) - 2*(sqrtS**2 + p1k2*x1) + 
     &                   muu2**2*(-4*p1k1 + sqrtS**2*(2*(2 + muu1**2 + muu3**2 + muu4**2 - x1) - x2)) + 2*p1k1*x2 + sqrtS**2*(2*(-1 + muu1**2)*(muu3**2 + muu4**2 - x1) + x2 - 3*muu1**2*x2))*
     &                 FR(indj,i,tt,gen)) + muu1*ELc(indi,i,tt,gen)*(muu2*sqrtS**2*
     &                 (muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + 
     &                   muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x1)*(-1 + x1 + x2))*FL(indj,i,tt,gen) + 
     &                muu4*(2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                   sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                      sqrtS**2*(-((1 + muu4**2 - x1)*(1 + x1)) + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*FR(indj,i,tt,gen))) + 
     &          dconjg(DR(j))*(2*muu1*muu3*ELc(indi,i,tt,gen)*(2*muu2*muu4*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*
     &                 FL(indj,i,tt,gen) + (-2*sqrtS**2*(p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) + sqrtS**2*(muu4**2*(-1 + muu1**2 + x1) + muu2**2*(-1 - muu3**2 + x1))) + 
     &                   4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**4*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1)*x2)*FR(indj,i,tt,gen)) + 
     &             ER(indi,i,tt,gen)*(4*muu2*muu4*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 + muu3**2 + muu4**2 - x1)) + 
     &                   muu1**2*sqrtS**4*(1 + muu4**2 - x2))*FL(indj,i,tt,gen) + 
     &                (2*p1k1*(8*p1k2 + sqrtS**2*((-3 + 2*muu2**2 + muu3**2 - x1)*x2 + muu4**2*(4 - 4*muu2**2 + x2))) + 
     &                   sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-1 - 4*muu1**2 + 2*muu2**2 - muu4**2 + x1) + 
     &                      sqrtS**2*(-((-1 + muu2**2)*(-1 + 2*muu2**2 - muu4**2 + x1)*(1 + muu4**2 - x1 - x2)) + 
     &                         muu1**2*(muu4**4 + (-1 + x1)**2 + 2*muu2**2*(1 + muu4**2 - x1 - x2) + (-1 + 5*x1)*x2 - muu4**2*(-2 + 6*x1 + x2))) + 
     &                      muu3**2*(p1k2*(8*muu1**2 - 2*(4 + x1)) + sqrtS**2*(muu1**2*(2 + 2*muu2**2 + 6*muu4**2 - 2*x1 - 5*x2) - (-1 + muu2**2)*(2*(-1 + muu2**2 - 3*muu4**2 + x1) + 5*x2)))
     &                      ))*FR(indj,i,tt,gen))))))/
     &  (2.*(1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmH(j))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen))))  )

        if (cutoff.lt.zerocheck) then

** TSchannel dipole from Higgs 312

       ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal( (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (((4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &           (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &              x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*
     &         (-(dconjg(CL(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 - 
     &                   (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) - 
     &                   4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                   (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*FL(indj,i,tt,gen)
     &                 + 4*muu2*muu4*sqrtS**6*(-4*muu4**2 + x2**2)*FR(indj,i,tt,gen))) + 
     &           dconjg(CR(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*(4*muu2*muu4*sqrtS**6*(4*muu4**2 - x2**2)*FL(indj,i,tt,gen) + 
     &              (4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 + 
     &                 (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) + 
     &                 4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*FR(indj,i,tt,gen)))
     &         )/(4.*sqrtS**8*(-4*muu4**2 + x2**2)**2) - ((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 - 
     &           (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) - 
     &           4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &           (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &         (muu1*muu4*sqrtS**2*(dconjg(CL(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)) + 
     &           ((-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                   x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*
     &              (dconjg(CL(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(4.*sqrtS**4*(4*muu4**2 - x2**2))
     &           ))/(sqrtS**4*(4*muu4**2 - x2**2)) + sqrtS**2*(-((muu2*muu3*
     &              (-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                   x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*
     &              (dconjg(CR(j))*dconjg(DL(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DR(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(sqrtS**4*(4*muu4**2 - x2**2)))
     &          + (muu1*muu3*(4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 + 
     &              (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) + 
     &              4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &              (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &            (dconjg(CL(j))*dconjg(DL(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(sqrtS**4*(4*muu4**2 - x2**2)) + 
     &         sqrtS**2*(dconjg(CL(j))*((-1 + muu1**2 + muu2**2)*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + 
     &               2*muu1*muu2*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen)) + 
     &            dconjg(CR(j))*(2*muu1*muu2*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + 
     &               (-1 + muu1**2 + muu2**2)*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*dconjg(DenSmH(j))*(-((muu1**2 + muu3**2)*sqrtS**2) + 
     &      (4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &         (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &            x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/
     &       (2.*sqrtS**4*(4*muu4**2 - x2**2)) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))) )

** TSchannel dipole from Higgs 321

       ChiChi2QQ_RadiationTS = ChiChi2QQ_RadiationTS - 2*dreal( (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) + (Sqrt((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*
     &         (-3*muu3**4 + muu4**4 + muu4**2*(x1 - x2) - (-1 + x1)*(-1 + 2*x1 + x2) + muu3**2*(-4 - 2*muu4**2 + 5*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (((((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &           (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &           ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &         (dconjg(CL(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*((4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &                 4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                 (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*FL(indj,i,tt,gen) + 4*muu2*muu4*sqrtS**6*(4*muu3**2 - x1**2)*FR(indj,i,tt,gen)) + 
     &           dconjg(CR(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*(4*muu2*muu4*sqrtS**6*(4*muu3**2 - x1**2)*FL(indj,i,tt,gen) + 
     &              (4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &                 4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                 (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*FR(indj,i,tt,gen))))/(sqrtS**4*(4*muu3**2 - x1**2)) - 
     &      ((4*(1 - muu1**2 + muu2**2)*muu3**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 + 
     &           4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*
     &         (muu1*muu4*sqrtS**2*(dconjg(CL(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)) - 
     &           ((-4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                   x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &              (dconjg(CL(j))*dconjg(DL(j))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(4.*sqrtS**4*(4*muu3**2 - x1**2))
     &           ))/(sqrtS**4*(4*muu3**2 - x1**2)) + sqrtS**2*(-((muu2*muu3*
     &              (4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                   x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &              (dconjg(CR(j))*dconjg(DL(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DR(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(sqrtS**4*(4*muu3**2 - x1**2)))
     &          + (muu1*muu3*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &              4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &              (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*
     &            (dconjg(CL(j))*dconjg(DL(j))*ER(indi,i,tt,gen)*FL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*ELc(indi,i,tt,gen)*FR(indj,i,tt,gen)))/(sqrtS**4*(4*muu3**2 - x1**2)) + 
     &         sqrtS**2*(dconjg(CL(j))*((-1 + muu1**2 + muu2**2)*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + 
     &               2*muu1*muu2*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen)) + 
     &            dconjg(CR(j))*(2*muu1*muu2*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*ELc(indi,i,tt,gen)*FL(indj,i,tt,gen) + 
     &               (-1 + muu1**2 + muu2**2)*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*ER(indi,i,tt,gen)*FR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dconjg(DenSmH(j))*(-((muu1**2 + muu3**2)*sqrtS**2) + 
     &      2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))) + (MSf2(i,itypeT,gen)-cI*MSf(i,itypeT,gen)*WSf(i,itypeT,gen)))) )

        endif

	  enddo

	enddo

	end

************************************************

	double precision function ChiChi2QQ_RadiationUS()

c#include "mykinematics.h"
c#include "model.h"
c#include "gencoupl.h"

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j,imin,imax,jmin,jmax
	complex*16 DenSmZ, DenSmH(6)
	double precision DenSmZ1, DenSmH1(6)

c       wZ = 0d0
c	wW = 0d0
	
c	wHk0(1) = 0d0
c	wHk0(2) = 0d0
c	wHk0(3) = 0d0
c	wHk0(4) = 0d0

c	wHk(1) = 0d0
c	wHk(2) = 0d0
c	wHk(3) = 0d0
c	wHk(4) = 0d0
c	wHk(5) = 0d0
c	wHk(6) = 0d0

	if (icha.eq.0) then
		itypeT = tt
		itypeU = tt
	else if (icha.eq.1) then
* tt is defined as normally as down-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt
		else if (tt.eq.4) then 
		      itypeT = tt
		      itypeU = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt+1
		else if (tt.eq.4) then 
		      itypeT = tt-1
		      itypeU = tt-1
		endif	
	endif	
    

        if (icha.eq.1) then
c only charged Higgs exchange		
	  imin = 5
	  imax = 6 
	  jmin = 5
	  jmax = 6
	else
c exchange of neutral Higgses 		
	  imin = 1
	  imax = 4 
	  jmin = 1
	  jmax = 4	  	
	endif

	if (icha.eq.1) then

	DenSmZ = (sqrtS**2 - MW2) + cI*MW*wW
        DenSmZc = (sqrtS**2 - MW2) - cI*MW*wW
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MW2)
		
	else
		
	DenSmZ = (sqrtS**2 - MZ2) + cI*MZ*wZ
        DenSmZc = (sqrtS**2 - MZ2) - cI*MZ*wZ
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MZ2)

	endif

	do i = imin, imax
	  DenSmH(i) = (sqrtS**2 - MHk0(i)**2) + cI*MHk0(i)*wHk0(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk0(i)**2)
	enddo

	ChiChi2QQ_RadiationUS = 0d0

	do i = 1, 2
** USchannel bremsstrahlung 11

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal( (4*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((4*p1k1*p1k2*(-1 + x2) + 4*p1k2**2*(-1 + 2*muu3**2 + x2)
     &            - 2*p1k2*sqrtS**2*(-1 + x2)*(-1 + 2*muu3**2 + x1 + x2) + 
     &            sqrtS**4*(muu3**2*(1 + 2*muu3**2 - x2) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*HL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(muu3**2*(8*p1k2 - 2*sqrtS**2*(-3 + x2)) + 2*p1k1*(-1 + x2) - sqrtS**2*(-1 + x2)*(-3 + 2*x1 + 3*x2))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(muu3**2*(8*p1k2 - 2*sqrtS**2*(-3 + x2)) + 2*p1k1*(-1 + x2) - sqrtS**2*(-1 + x2)*(-3 + 2*x1 + 3*x2))*HL(indj,i,tt,gen) + 
     &         (4*p1k1*p1k2*(-1 + x2) + 4*p1k2**2*(-1 + 2*muu3**2 + x2) - 2*p1k2*sqrtS**2*(-1 + x2)*(-1 + 2*muu3**2 + x1 + x2) + 
     &            sqrtS**4*(muu3**2*(1 + 2*muu3**2 - x2) + muu1**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*HR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x2)**2*(-2*p1k2 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))  )

	endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (4*CF*gs**2*NC*(dconjg(AL)*(-2*dconjg(BL)*(p1k2*GL(indi,i,tt,gen) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*
     &          (-2*muu2*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HL(indj,i,tt,gen) + 2*p1k1*HR(indj,i,tt,gen) + 
     &            (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 - 2*(p1k1 + p1k2)*x2 + sqrtS**2*(x1*(-1 + x2) + x2*(-2 + muu1**2 - muu2**2 + x2)) + 
     &               muu4**2*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2)) + muu3**2*(-2*p1k1 - 6*p1k2 + sqrtS**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*
     &             HR(indj,i,tt,gen)) + sqrtS**2*dconjg(BR)*(muu4*GL(indi,i,tt,gen)*
     &             (-(muu2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*HL(indj,i,tt,gen)) - 
     &               (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*HL(indj,i,tt,gen) + 
     &               muu3*(1 + muu3**2 + muu4**2 - x2)*(-2*p1k2 + sqrtS**2*x2)*HR(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(2*dconjg(BR)*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen))*
     &          (((-1 + muu1**2 - muu2**2)*sqrtS**2 + 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 2*p1k2*(-1 + 3*muu3**2 - muu4**2 + x2) - 
     &               sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*HL(indj,i,tt,gen) + 
     &            2*muu2*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu1*GL(indi,i,tt,gen)*(muu3*(1 + muu3**2 + muu4**2 - x2)*(-2*p1k2 + sqrtS**2*x2)*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HL(indj,i,tt,gen)) - 
     &               muu2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)**2*dreal(DenSmZ)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 12

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   (2*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((-4*p1k2**2*x1 + sqrtS**4*(2*(muu1**2 + muu3**2) + x1)
     &            *(1 + 2*muu3**2 - x1 - x2) - 4*p1k1**2*x2 + 
     &            2*p1k2*sqrtS**2*(1 - x2 + x1*(-3 + 2*muu3**2 + x1 + 2*x2)) + 2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(-1 + x1 + 2*muu3**2*(-2 + x2) + x2 + x1*x2)))*
     &          HL(indj,i,tt,gen) + 2*muu1*muu3*sqrtS**2*(-2*p1k2*x1 + 2*p1k1*(1 - 4*muu3**2 - x1 + x2) + sqrtS**2*(4 + 2*muu3**2*(3 + x1) - 4*x2 + x1*(-6 + x1 + 2*x2)))*HR(indj,i,tt,gen))
     &        + GL(indi,i,tt,gen)*(2*muu1*muu3*sqrtS**2*(-2*p1k2*x1 + 2*p1k1*(1 - 4*muu3**2 - x1 + x2) + sqrtS**2*(4 + 2*muu3**2*(3 + x1) - 4*x2 + x1*(-6 + x1 + 2*x2)))*
     &          HL(indj,i,tt,gen) + (-4*p1k2**2*x1 + sqrtS**4*(2*(muu1**2 + muu3**2) + x1)*(1 + 2*muu3**2 - x1 - x2) - 4*p1k1**2*x2 + 
     &            2*p1k2*sqrtS**2*(1 - x2 + x1*(-3 + 2*muu3**2 + x1 + 2*x2)) + 2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(-1 + x1 + 2*muu3**2*(-2 + x2) + x2 + x1*x2)))*
     &          HR(indj,i,tt,gen))))/(sqrtS**4*(-1 + x1)*(-1 + x2)*(-2*p1k2 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))  )

	endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (2*CF*gs**2*NC*(dconjg(AR)*(dconjg(BR)*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             ((2*(-1 + muu3**2 + muu4**2)*p1k1 + (-1 + muu1**2 - muu2**2)*sqrtS**2 + 2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) + 
     &                  sqrtS**2*(muu4**2*(-1 + muu1**2 - muu2**2 - x1) + 2*x1 + muu2**2*x1 + x2 + muu2**2*x2 - x1*x2 - muu1**2*(x1 + x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HR(indj,i,tt,gen)) - 
     &            GR(indi,i,tt,gen)*((8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 
     &                  2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + muu4**2*
     &                   (8*p1k1**2 + 2*(1 - muu1**2 + muu2**2)*p1k2*sqrtS**2 + (1 + muu1**2 - muu2**2)*sqrtS**4*x1 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1)) - 4*p1k1**2*x2 - 
     &                  sqrtS**2*(2*p1k2*(1 - muu1**2 + muu2**2 + 2*muu3**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) - (1 + muu1**2 - muu2**2 + x1)*x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(2*(-1 + muu3**2)*p1k2 + muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - x1 - x2) + 
     &                  2*p1k1*x2)*HR(indj,i,tt,gen))) + dconjg(BL)*(muu1*GL(indi,i,tt,gen)*
     &             (muu3*sqrtS**2*(2*p1k2*(-2 + 4*muu4**2 + x1) + 2*p1k1*(-2*muu3**2 + 2*(-1 + muu4**2 + x1) + x2) + 
     &                  sqrtS**2*((1 + muu3**2 - muu4**2 - x1)*(-3 + 3*muu1**2 - 3*muu2**2 + 2*x1) + (3 - muu1**2 + muu2**2 - 4*muu4**2 - 2*x1)*x2))*HL(indj,i,tt,gen) + 
     &               2*muu2*sqrtS**4*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*(2*muu3*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*
     &                   (muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                  2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**2*(8*p1k2 - 3*(1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 - 2*p1k1*(-2 + x2) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu4**2 + x1 + x2))*
     &                HR(indj,i,tt,gen)))) + dconjg(AL)*(dconjg(BR)*(muu1*sqrtS**2*GR(indi,i,tt,gen)*
     &             (2*muu2*sqrtS**2*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*HL(indj,i,tt,gen) + 
     &               muu3*(2*p1k2*(-2 + 4*muu4**2 + x1) + 2*p1k1*(-2*muu3**2 + 2*(-1 + muu4**2 + x1) + x2) + 
     &                  sqrtS**2*((1 + muu3**2 - muu4**2 - x1)*(-3 + 3*muu1**2 - 3*muu2**2 + 2*x1) + (3 - muu1**2 + muu2**2 - 4*muu4**2 - 2*x1)*x2))*HR(indj,i,tt,gen)) + 
     &            muu4*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**2*(8*p1k2 - 3*(1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 - 2*p1k1*(-2 + x2) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu4**2 + x1 + x2))*HL(indj,i,tt,gen) + 
     &               2*muu3*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*(muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                  2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*HR(indj,i,tt,gen))) + 
     &         dconjg(BL)*(2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HL(indj,i,tt,gen) - 2*p1k1*HR(indj,i,tt,gen) - 
     &               (-2*p1k2*x1 + muu4**2*(-2*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1)) + muu3**2*(-2*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) + 
     &                  sqrtS**2*(1 - 2*x1 + (-1 + x1)*x2 + muu1**2*(-1 + x1 + x2) - muu2**2*(-1 + x1 + x2)))*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(-2*muu2*muu3*sqrtS**2*(2*(-1 + muu3**2)*p1k2 + muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - x1 - x2) + 2*p1k1*x2)*HL(indj,i,tt,gen) - 
     &               (8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 
     &                  2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + muu4**2*
     &                   (8*p1k1**2 + 2*(1 - muu1**2 + muu2**2)*p1k2*sqrtS**2 + (1 + muu1**2 - muu2**2)*sqrtS**4*x1 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1)) - 4*p1k1**2*x2 - 
     &                  sqrtS**2*(2*p1k2*(1 - muu1**2 + muu2**2 + 2*muu3**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) - (1 + muu1**2 - muu2**2 + x1)*x2)))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmZ)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 21

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   (-2*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((-4*p1k2**2*x1 - 4*p1k1**2*x2 + 2*p1k2*sqrtS**2*(1 + 2*muu3**2*(2 + x1)
     &            - x2 + x1*(-3 + x1 + 2*x2)) + 
     &            2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(-1 + x1 + x2 + 2*muu3**2*x2 + x1*x2)) + 
     &            sqrtS**4*(4*muu3**4 + x1 - x1*(x1 + x2) - 2*muu3**2*(-1 + x1 + 2*x2) + muu1**2*(4*muu3**2 - 2*(-1 + x1 + x2))))*HL(indj,i,tt,gen) + 
     &         2*muu1*muu3*sqrtS**2*(sqrtS**2*(4 - 2*muu3**2*(-3 + x2) + 2*x1*(-2 + x2) - 5*x2) + 2*p1k1*x2 + 2*p1k2*(-1 + 4*muu3**2 - x1 + x2))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(2*muu1*muu3*sqrtS**2*(sqrtS**2*(4 - 2*muu3**2*(-3 + x2) + 2*x1*(-2 + x2) - 5*x2) + 2*p1k1*x2 + 2*p1k2*(-1 + 4*muu3**2 - x1 + x2))*HL(indj,i,tt,gen) + 
     &         (-4*p1k2**2*x1 - 4*p1k1**2*x2 + 2*p1k2*sqrtS**2*(1 + 2*muu3**2*(2 + x1) - x2 + x1*(-3 + x1 + 2*x2)) + 
     &            2*p1k1*(-2*p1k2*(-2 + 4*muu3**2 + x1 + x2) + sqrtS**2*(-1 + x1 + x2 + 2*muu3**2*x2 + x1*x2)) + 
     &            sqrtS**4*(4*muu3**4 + x1 - x1*(x1 + x2) - 2*muu3**2*(-1 + x1 + 2*x2) + muu1**2*(4*muu3**2 - 2*(-1 + x1 + x2))))*HR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)*(-1 + x2)*(-2*p1k1 - sqrtS**2*(muu1**2 + muu3**2 - x1) + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))  )

	endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (2*CF*gs**2*NC*(dconjg(AR)*(dconjg(BR)*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             ((2*(-1 + muu3**2 + muu4**2)*p1k1 + (-1 + muu1**2 - muu2**2)*sqrtS**2 + 2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) + 
     &                  sqrtS**2*(muu4**2*(-1 + muu1**2 - muu2**2 - x1) + 2*x1 + muu2**2*x1 + x2 + muu2**2*x2 - x1*x2 - muu1**2*(x1 + x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HR(indj,i,tt,gen)) - 
     &            GR(indi,i,tt,gen)*((8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 
     &                  2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + muu4**2*
     &                   (8*p1k1**2 + 2*(1 - muu1**2 + muu2**2)*p1k2*sqrtS**2 + (1 + muu1**2 - muu2**2)*sqrtS**4*x1 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1)) - 4*p1k1**2*x2 - 
     &                  sqrtS**2*(2*p1k2*(1 - muu1**2 + muu2**2 + 2*muu3**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) - (1 + muu1**2 - muu2**2 + x1)*x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(2*(-1 + muu3**2)*p1k2 + muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - x1 - x2) + 
     &                  2*p1k1*x2)*HR(indj,i,tt,gen))) + dconjg(BL)*(muu1*GL(indi,i,tt,gen)*
     &             (muu3*sqrtS**2*(2*p1k2*(-2 + x1) + muu4**2*(-8*p1k1 + sqrtS**2*(-3 + 3*muu1**2 - 3*muu2**2 + 4*x1)) + 2*p1k1*x2 + 
     &                  sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2))*HL(indj,i,tt,gen) + 
     &               2*muu2*sqrtS**4*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*(2*muu3*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 
     &                  2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2))*HR(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(dconjg(BR)*(muu1*sqrtS**2*GR(indi,i,tt,gen)*(2*muu2*sqrtS**2*(-1 + x1 - muu4**2*(1 + x1) + muu3**2*(-1 + 4*muu4**2 - x2) + x2)*HL(indj,i,tt,gen) + 
     &               muu3*(2*p1k2*(-2 + x1) + muu4**2*(-8*p1k1 + sqrtS**2*(-3 + 3*muu1**2 - 3*muu2**2 + 4*x1)) + 2*p1k1*x2 + 
     &                  sqrtS**2*((-1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) + (3 - muu1**2 + muu2**2 - 2*x1)*x2))*HR(indj,i,tt,gen)) + 
     &            muu4*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2))*HL(indj,i,tt,gen) + 
     &               2*muu3*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                  sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*HR(indj,i,tt,gen))) + 
     &         dconjg(BL)*(2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HL(indj,i,tt,gen) - 2*p1k1*HR(indj,i,tt,gen) - 
     &               (-2*p1k2*x1 + muu4**2*(-2*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1)) + muu3**2*(-2*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + x1 - 2*x2)) + 
     &                  sqrtS**2*(1 - 2*x1 + (-1 + x1)*x2 + muu1**2*(-1 + x1 + x2) - muu2**2*(-1 + x1 + x2)))*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(-2*muu2*muu3*sqrtS**2*(2*(-1 + muu3**2)*p1k2 + muu4**2*(-4*p1k1 + 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - x1 - x2) + 2*p1k1*x2)*HL(indj,i,tt,gen) - 
     &               (8*muu3**2*p1k2**2 + 2*p1k2*sqrtS**2 - 2*muu1**2*p1k2*sqrtS**2 + 2*muu2**2*p1k2*sqrtS**2 + 2*muu3**2*p1k2*sqrtS**2 - 2*muu1**2*muu3**2*p1k2*sqrtS**2 + 
     &                  2*muu2**2*muu3**2*p1k2*sqrtS**2 - 4*p1k2**2*x1 - 6*p1k2*sqrtS**2*x1 + 2*muu1**2*p1k2*sqrtS**2*x1 - 2*muu2**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 + muu1**2*sqrtS**4*x1 - 
     &                  muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + muu1**2*muu3**2*sqrtS**4*x1 - muu2**2*muu3**2*sqrtS**4*x1 + 2*p1k2*sqrtS**2*x1**2 - sqrtS**4*x1**2 - 
     &                  muu1**2*sqrtS**4*x1**2 + muu2**2*sqrtS**4*x1**2 + muu4**2*
     &                   (8*p1k1**2 + 2*(1 - muu1**2 + muu2**2)*p1k2*sqrtS**2 + (1 + muu1**2 - muu2**2)*sqrtS**4*x1 - 2*p1k1*sqrtS**2*(1 + muu1**2 - muu2**2 + 2*x1)) - 4*p1k1**2*x2 - 
     &                  sqrtS**2*(2*p1k2*(1 - muu1**2 + muu2**2 + 2*muu3**2 - 2*x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*x1)*x2 - 
     &                  2*p1k1*(2*p1k2*(-2 + x1 + x2) + sqrtS**2*((1 + muu1**2 - muu2**2)*(1 + muu3**2 - x1) - (1 + muu1**2 - muu2**2 + x1)*x2)))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmZ)*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 22

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   (4*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((4*p1k1**2*(-1 + 2*muu3**2 + x1) + 2*p1k1*(2*p1k2*(-1 + x1) - sqrtS**2*(1 + x1)
     &            *(-1 + 2*muu3**2 + x1)) + 
     &            sqrtS**2*(2*muu3**4*sqrtS**2 - (2*p1k2 - sqrtS**2)*(-1 + x1)*x1 + muu3**2*sqrtS**2*(1 + x1) + muu1**2*sqrtS**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*HL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(muu3**2*(-8*p1k1 + 2*sqrtS**2*(3 + x1)) - (-1 + x1)*(2*p1k2 + sqrtS**2*(-3 + 3*x1 + x2)))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(muu3**2*(-8*p1k1 + 2*sqrtS**2*(3 + x1)) - (-1 + x1)*(2*p1k2 + sqrtS**2*(-3 + 3*x1 + x2)))*HL(indj,i,tt,gen) + 
     &         (4*p1k1**2*(-1 + 2*muu3**2 + x1) + 2*p1k1*(2*p1k2*(-1 + x1) - sqrtS**2*(1 + x1)*(-1 + 2*muu3**2 + x1)) + 
     &            sqrtS**2*(2*muu3**4*sqrtS**2 - (2*p1k2 - sqrtS**2)*(-1 + x1)*x1 + muu3**2*sqrtS**2*(1 + x1) + muu1**2*sqrtS**2*(-1 + 2*muu3**2 + x1 + x2 - x1*x2)))*HR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)**2*(2*p1k1 + sqrtS**2*(muu1**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

	endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (4*CF*gs**2*NC*(dconjg(AR)*(dconjg(BR)*(2*muu1*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*GL(indi,i,tt,gen) + 2*p1k1*GR(indi,i,tt,gen) + 
     &            (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2 + muu4**2*(-6*p1k1 - 2*p1k2 + 3*(1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &               muu3**2*(2*(p1k1 + p1k2) + (-1 - muu1**2 + muu2**2)*sqrtS**2) + (-2*(p1k1 + p1k2) + (1 + muu1**2 - muu2**2)*sqrtS**2)*x1)*GR(indi,i,tt,gen))*
     &          ((-2*p1k1 + sqrtS**2*x1)*HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu4*(1 + muu3**2 + muu4**2 - x1)*GR(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HL(indj,i,tt,gen)) + 2*muu2*p1k1*HR(indj,i,tt,gen)) + 
     &            muu1*GL(indi,i,tt,gen)*(-(muu3*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                    sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*HL(indj,i,tt,gen)) + 
     &               muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*HR(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(dconjg(BL)*((2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 - muu3**2 + 3*muu4**2 + x1))*
     &             GL(indi,i,tt,gen) + 2*muu1*muu4*sqrtS**2*(1 + muu3**2 + muu4**2 - x1)*GR(indi,i,tt,gen))*
     &          (2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + (-2*p1k1 + sqrtS**2*x1)*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu4*(1 + muu3**2 + muu4**2 - x1)*GL(indi,i,tt,gen)*(2*muu2*p1k1*HL(indj,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*HL(indj,i,tt,gen) - 
     &               muu3*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                  sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)**2*dreal(DenSmZ)*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 31

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal( (2*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((-2*p1k1*(8*p1k2**2 + (muu1**2 + muu3**2)*sqrtS**4*x2 - 2*p1k2*sqrtS**2*(-1 + 2*muu3**2 + x2)) + 
     &            sqrtS**2*(-4*muu3**4*sqrtS**4 + 2*p1k2*x1*(4*p1k2 - sqrtS**2*(-1 + x2)) + muu1**2*sqrtS**2*(sqrtS**2*x2*(-1 + x1 + x2) + 2*p1k2*(-2 + x1 + 2*x2)) + 
     &               muu3**2*(-8*p1k2**2 - 2*p1k2*sqrtS**2*(2 + x1 - 2*x2) + sqrtS**4*(-2 + x1 + 2*x2 - 4*muu1**2*(-1 + x1 + x2)))))*HL(indj,i,tt,gen) + 
     &         muu1*muu3*sqrtS**2*(8*p1k2**2 + p1k1*(-8*p1k2 + 2*sqrtS**2*(-3 + x2)) + 2*p1k2*sqrtS**2*(-1 - 8*muu3**2 + 4*x1 + x2) + 
     &            sqrtS**4*(-2 - 2*muu3**2*(4 + x1 - x2) + x2 + 2*muu1**2*x2 + x2**2 - 2*x1*(-2 + muu1**2 + x2)))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(8*p1k2**2 + p1k1*(-8*p1k2 + 2*sqrtS**2*(-3 + x2)) + 2*p1k2*sqrtS**2*(-1 - 8*muu3**2 + 4*x1 + x2) + 
     &            sqrtS**4*(-2 - 2*muu3**2*(4 + x1 - x2) + x2 + 2*muu1**2*x2 + x2**2 - 2*x1*(-2 + muu1**2 + x2)))*HL(indj,i,tt,gen) + 
     &         (-2*p1k1*(8*p1k2**2 + (muu1**2 + muu3**2)*sqrtS**4*x2 - 2*p1k2*sqrtS**2*(-1 + 2*muu3**2 + x2)) + 
     &            sqrtS**2*(-4*muu3**4*sqrtS**4 + 2*p1k2*x1*(4*p1k2 - sqrtS**2*(-1 + x2)) + muu1**2*sqrtS**2*(sqrtS**2*x2*(-1 + x1 + x2) + 2*p1k2*(-2 + x1 + 2*x2)) + 
     &               muu3**2*(-8*p1k2**2 - 2*p1k2*sqrtS**2*(2 + x1 - 2*x2) + sqrtS**4*(-2 + x1 + 2*x2 - 4*muu1**2*(-1 + x1 + x2)))))*HR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x2)*(-2*p1k2 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(2*p1k1 + sqrtS**2*(muu1**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))  )

	endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (2*CF*gs**2*NC*(dconjg(AR)*(2*dconjg(BR)*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen))*
     &          (-((-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &                 sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*
     &               HL(indj,i,tt,gen)) + 2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu4*GR(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*(2*muu2**2*p1k2 + sqrtS**2 - p1k2*x1) + 
     &                  muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) - 2*p1k1*x2 + 
     &                  sqrtS**2*(-2*muu2**4 - 2*muu3**2 - 2*muu4**2 + x1 + muu2**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2) + 2*x2))*HL(indj,i,tt,gen) + 
     &               muu2*(8*p1k1*p1k2 - 2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k1*sqrtS**2 + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                HR(indj,i,tt,gen)) + muu1*GL(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                  sqrtS**2*(2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - (muu3**2 + muu4**2)*x2 + x2**2 + 
     &                        muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 + muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + 
     &                  muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x2)*(-1 + x1 + x2))*HR(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(2*dconjg(BL)*(p1k2*GL(indi,i,tt,gen) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*
     &          (2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*HL(indj,i,tt,gen) - 
     &            (-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &               sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*
     &             HR(indj,i,tt,gen)) + sqrtS**2*dconjg(BR)*(muu4*GL(indi,i,tt,gen)*
     &             (muu2*(8*p1k1*p1k2 - 2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k1*sqrtS**2 + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                HL(indj,i,tt,gen) + muu3*sqrtS**2*(-2*(2*muu2**2*p1k2 + sqrtS**2 - p1k2*x1) + muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) - 
     &                  2*p1k1*x2 + sqrtS**2*(-2*muu2**4 - 2*muu3**2 - 2*muu4**2 + x1 + muu2**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2) + 2*x2))*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 + muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + 
     &                  muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x2)*(-1 + x1 + x2))*HL(indj,i,tt,gen) + 
     &               muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + sqrtS**2*
     &                   (2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + sqrtS**2*
     &                      (-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - (muu3**2 + muu4**2)*x2 + x2**2 + 
     &                        muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmZ)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 32

        if (icha.eq.2) then
		
** US-channel photon (for cha cha -> b bbar)
** relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (-2*ApL*BpL*CF*gs**2*NC*(GR(indi,i,tt,gen)*((8*p1k1**2*(2*p1k2 - muu3**2*sqrtS**2) + 
     &            2*p1k1*sqrtS**2*(p1k2*(-2 + 4*muu3**2 - 6*x1) - sqrtS**2*(muu3**2*(-2 - 2*x1 + x2) + muu1**2*(-2 + 2*x1 + x2))) + 
     &            sqrtS**4*(-4*muu3**4*sqrtS**2 + x1*(2*p1k2*(1 + muu1**2 + x1) + muu1**2*sqrtS**2*(-3 + 3*x1 + x2)) + muu3**2*(-2*p1k2*x1 + sqrtS**2*(-2 + x2 - 4*muu1**2*(-1 + x1 + x2))))
     &            )*HL(indj,i,tt,gen) + muu1*muu3*sqrtS**2*(8*p1k1**2 + 2*p1k2*sqrtS**2*(3 + x1) + 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 8*muu3**2 - 5*x1 - 2*x2)) + 
     &            sqrtS**4*(-2 + 2*muu1**2*(x1 - x2) + x2 - 2*muu3**2*(4 + 3*x1 + x2) + x1*(4*x1 + x2)))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(muu1*muu3*sqrtS**2*(8*p1k1**2 + 2*p1k2*sqrtS**2*(3 + x1) + 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 8*muu3**2 - 5*x1 - 2*x2)) + 
     &            sqrtS**4*(-2 + 2*muu1**2*(x1 - x2) + x2 - 2*muu3**2*(4 + 3*x1 + x2) + x1*(4*x1 + x2)))*HL(indj,i,tt,gen) + 
     &         (8*p1k1**2*(2*p1k2 - muu3**2*sqrtS**2) + 2*p1k1*sqrtS**2*(p1k2*(-2 + 4*muu3**2 - 6*x1) - sqrtS**2*(muu3**2*(-2 - 2*x1 + x2) + muu1**2*(-2 + 2*x1 + x2))) + 
     &            sqrtS**4*(-4*muu3**4*sqrtS**2 + x1*(2*p1k2*(1 + muu1**2 + x1) + muu1**2*sqrtS**2*(-3 + 3*x1 + x2)) + muu3**2*(-2*p1k2*x1 + sqrtS**2*(-2 + x2 - 4*muu1**2*(-1 + x1 + x2))))
     &            )*HR(indj,i,tt,gen))))/
     &  (sqrtS**4*(-1 + x1)*(-2*p1k2 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(-2*p1k1 - sqrtS**2*(muu1**2 + muu3**2 - x1) + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

        endif	

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (-2*CF*gs**2*NC*(dconjg(AR)*(-2*dconjg(BR)*(muu1*muu4*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2
     &            *(-1 + muu1**2 - muu2**2 - 2*muu3**2 - 2*muu4**2 + 2*x1))*GL(indi,i,tt,gen) - 
     &            (p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + 
     &                  p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*GR(indi,i,tt,gen))*((2*p1k1 - sqrtS**2*x1)*HL(indj,i,tt,gen) - 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu4*GR(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*muu1**4*sqrtS**2 + 2*p1k2*x1 + muu1**2*(4*p1k1 + sqrtS**2*(2*(2 + muu2**2 + muu3**2 + muu4**2 - 2*x1) - x2)) + 
     &                  muu2**2*(-4*p1k2 + sqrtS**2*(2*(1 + muu3**2 + muu4**2 - x1) - x2)) - 2*p1k1*x2 + sqrtS**2*(-2*(1 + muu3**2 + muu4**2 - x1) + x2))*HL(indj,i,tt,gen) + 
     &               muu2*(sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*
     &                HR(indj,i,tt,gen)) + muu1*GL(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                  sqrtS**2*(2*p1k2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - (muu3**2 + muu4**2 - x1)*x2 + 
     &                        muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x1)*(-1 + x1 + x2) + 
     &                  muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HR(indj,i,tt,gen))))
     &        + dconjg(AL)*(-2*dconjg(BL)*((p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + 
     &               sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*GL(indi,i,tt,gen) + 
     &            muu1*muu4*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*GR(indi,i,tt,gen))*
     &          (2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + (-2*p1k1 + sqrtS**2*x1)*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BR)*(muu4*GL(indi,i,tt,gen)*(muu2*(sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*
     &                HL(indj,i,tt,gen) + muu3*sqrtS**2*(-2*muu1**4*sqrtS**2 + 2*p1k2*x1 + muu1**2*(4*p1k1 + sqrtS**2*(2*(2 + muu2**2 + muu3**2 + muu4**2 - 2*x1) - x2)) + 
     &                  muu2**2*(-4*p1k2 + sqrtS**2*(2*(1 + muu3**2 + muu4**2 - x1) - x2)) - 2*p1k1*x2 + sqrtS**2*(-2*(1 + muu3**2 + muu4**2 - x1) + x2))*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x1)*(-1 + x1 + x2) + 
     &                  muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HL(indj,i,tt,gen) + 
     &               muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + sqrtS**2*
     &                   (2*p1k2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - (muu3**2 + muu4**2 - x1)*x2 + 
     &                        muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dreal(DenSmZ)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

        if (cutoff.lt.zerocheck) then

** USchannel dipole Z 312

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal(   (4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (dconjg(AL)*(-4*dconjg(BL)*(((-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                    x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*GL(indi,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2)) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*
     &          (muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + ((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 - 
     &                 (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) - 
     &                 4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HR(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2))) + sqrtS**2*dconjg(BR)*
     &          (muu4*GL(indi,i,tt,gen)*((muu2*(-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                       x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HL(indj,i,tt,gen)
     &                  )/(2.*sqrtS**4*(4*muu4**2 - x2**2)) + (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*HL(indj,i,tt,gen) - 
     &               (muu3*(4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 + 
     &                    (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) + 
     &                    4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                    (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HR(indj,i,tt,gen)
     &                  )/(2.*sqrtS**4*(4*muu4**2 - x2**2))))) + dconjg(AR)*
     &       (-4*dconjg(BR)*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + ((-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                    x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*GR(indi,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2)))*(((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 - 
     &                 (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) - 
     &                 4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HL(indj,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu4**2 - x2**2)) + muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu1*GL(indi,i,tt,gen)*(-(muu3*(4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 + 
     &                     (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) + 
     &                     4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                     (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &                   HL(indj,i,tt,gen))/(2.*sqrtS**4*(4*muu4**2 - x2**2)) + muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HL(indj,i,tt,gen) + 
     &               (muu2*(-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                       x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HR(indj,i,tt,gen)
     &                  )/(2.*sqrtS**4*(4*muu4**2 - x2**2)))))))/
     &  (sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*dreal(DenSmZ)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 + (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) + 
     &         4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))/
     &       (2.*sqrtS**4*(4*muu4**2 - x2**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel dipole Z 321

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal(   (4*CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) + (Sqrt((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*
     &         (-3*muu3**4 + muu4**4 + muu4**2*(x1 - x2) - (-1 + x1)*(-1 + 2*x1 + x2) + muu3**2*(-4 - 2*muu4**2 + 5*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (dconjg(AL)*(-4*dconjg(BL)*(((4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                    x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*GL(indi,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu3**2 - x1**2)) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*
     &          (muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) - ((4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &                 4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                 (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*HR(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu3**2 - x1**2))) + 
     &         sqrtS**2*dconjg(BR)*(muu4*GL(indi,i,tt,gen)*((muu2*(4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                       x1*((-1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &                  HL(indj,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) + (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*HL(indj,i,tt,gen) - 
     &               (muu3*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &                    4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                    (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*HR(indj,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2))))) + 
     &      dconjg(AR)*(-4*dconjg(BR)*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + 
     &            ((4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                 (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                    x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*GR(indi,i,tt,gen))/
     &             (4.*sqrtS**4*(4*muu3**2 - x1**2)))*(-((4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &                  4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                  (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*HL(indj,i,tt,gen))/(4.*sqrtS**4*(4*muu3**2 - x1**2)) + muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &         sqrtS**2*dconjg(BL)*(muu1*GL(indi,i,tt,gen)*(-(muu3*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &                     4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                     (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                     (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*HL(indj,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) + 
     &               muu2*(-1 + muu3**2 + muu4**2)*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HL(indj,i,tt,gen) + 
     &               (muu2*(4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(-4 + 4*muu4**2 + x1**2) + 
     &                       x1*((-1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &                  HR(indj,i,tt,gen))/(2.*sqrtS**4*(4*muu3**2 - x1**2)))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dreal(DenSmZ)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &         4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &         (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )




       if(icha.eq.2) then

** USchannel dipole with G 312

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal(   (4*ApL*BpL*CF*GS**2*NC*(-2*muu3**2*(4*Sqrt(1 - 4*muu3**2) - 8*Sqrt(1 - 4*muu3**2)*x2 + 3*Sqrt(1 - 4*muu3**2)*x2**2 
     &         + Sqrt(1 - 4*muu3**2)*x1*(-2 + 3*x2) - 
     &         4*Sqrt(-4*muu3**2 + x2**2) + 4*x2*Sqrt(-4*muu3**2 + x2**2)) + 
     &      (-1 + x2)*(Sqrt(1 - 4*muu3**2)*(x1**2 + 3*x1*(-1 + x2) + x2*(-5 + 2*x2)) + 2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x2**2))) + 
     &      4*muu3**4*(Sqrt(1 - 4*muu3**2)*x1 - 2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x2**2)) + x2*(Sqrt(1 - 4*muu3**2) + 2*Sqrt(-4*muu3**2 + x2**2))))*
     &    (GR(indi,i,tt,gen)*((8*(-1 + 4*muu3**2)*p1k2**2*sqrtS**2 + 4*(1 - 4*muu3**2)*p1k2*sqrtS**4*x2 - 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + sqrtS**6*((2 + 8*muu1**2)*muu3**2 + 8*muu3**4 - (1 + 2*muu1**2)*x2**2))*HL(indj,i,tt,gen)
     &           + 2*muu1*muu3*(12*muu3**2*sqrtS**6 - 3*sqrtS**6*x2**2 - 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)))*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(2*muu1*muu3*(12*muu3**2*sqrtS**6 - 3*sqrtS**6*x2**2 - 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)))*HL(indj,i,tt,gen) + 
     &         (8*(-1 + 4*muu3**2)*p1k2**2*sqrtS**2 + 4*(1 - 4*muu3**2)*p1k2*sqrtS**4*x2 - 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + sqrtS**6*((2 + 8*muu1**2)*muu3**2 + 8*muu3**4 - (1 + 2*muu1**2)*x2**2))*HR(indj,i,tt,gen)
     &         )))/((-1 + 2*muu3**2)*sqrtS**2*(-1 + x2)**2*(-2 + x1 + x2)*Sqrt(-4*muu3**2 + x2**2)*
     &    (8*muu3**4*sqrtS**6 + (1 - 2*muu1**2)*sqrtS**6*x2**2 + 2*muu3**2*sqrtS**6*(-2 + 4*muu1**2 - x2**2) + 4*p1k2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) - 
     &      sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu3**2 + x2**2)) + 2*sqrtS**4*(-4*muu3**2 + x2**2)*(MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )


** USchannel dipole with G 321

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal( (4*ApL*BpL*CF*GS**2*NC*(4*muu3**4*(-2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x1**2)) + x1*(Sqrt(1 - 4*muu3**2) 
     &      + 2*Sqrt(-4*muu3**2 + x1**2)) + Sqrt(1 - 4*muu3**2)*x2) - 
     &      2*muu3**2*(4*Sqrt(1 - 4*muu3**2) + 3*Sqrt(1 - 4*muu3**2)*x1**2 - 4*Sqrt(-4*muu3**2 + x1**2) - 2*Sqrt(1 - 4*muu3**2)*x2 + 
     &         x1*(-8*Sqrt(1 - 4*muu3**2) + 4*Sqrt(-4*muu3**2 + x1**2) + 3*Sqrt(1 - 4*muu3**2)*x2)) + 
     &      (-1 + x1)*(2*(Sqrt(1 - 4*muu3**2) + Sqrt(-4*muu3**2 + x1**2)) + Sqrt(1 - 4*muu3**2)*(x1*(-5 + 2*x1) + 3*(-1 + x1)*x2 + x2**2)))*
     &    (GL(indi,i,tt,gen)*(2*muu1*muu3*(12*muu3**2*sqrtS**6 + 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &            sqrtS**2*x1*(3*sqrtS**4*x1 + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*HL(indj,i,tt,gen) + 
     &         (8*(-1 + 4*muu3**2)*p1k1**2*sqrtS**2 + 4*(1 - 4*muu3**2)*p1k1*sqrtS**4*x1 + 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &            sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + sqrtS**6*((2 + 8*muu1**2)*muu3**2 + 8*muu3**4 - (1 + 2*muu1**2)*x1**2))*HR(indj,i,tt,gen)
     &         ) + GR(indi,i,tt,gen)*((8*(-1 + 4*muu3**2)*p1k1**2*sqrtS**2 + 4*(1 - 4*muu3**2)*p1k1*sqrtS**4*x1 + 
     &            4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &            sqrtS**6*((2 + 8*muu1**2)*muu3**2 + 8*muu3**4 - (1 + 2*muu1**2)*x1**2))*HL(indj,i,tt,gen) + 
     &         2*muu1*muu3*(12*muu3**2*sqrtS**6 + 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &            sqrtS**2*x1*(3*sqrtS**4*x1 + Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*HR(indj,i,tt,gen))))/
     &  ((-1 + 2*muu3**2)*sqrtS**2*(-1 + x1)**2*Sqrt(-4*muu3**2 + x1**2)*(-2 + x1 + x2)*
     &    (8*muu3**4*sqrtS**6 + (1 - 2*muu1**2)*sqrtS**6*x1**2 + 2*muu3**2*sqrtS**6*(-2 + 4*muu1**2 - x1**2) - 4*p1k1*Sqrt((1 - 4*muu3**2)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &      sqrtS**2*Sqrt((1 - 4*muu3**2)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 2*sqrtS**4*(-4*muu3**2 + x1**2)*(MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )
 
        endif

        endif


	  do j = jmin, jmax

** USchannel bremsstrahlung 13

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  -((CF*gs**2*NC*(dconjg(CR(j))*(dconjg(DR(j))*(-4*muu2*muu3*(1 + muu3**2 + muu4**2 - x2)
     &              *(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen))*HL(indj,i,tt,gen) + 
     &              (2*muu1*muu4*((1 - muu1**2 + muu2**2)*sqrtS**2 + 2*p1k2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*GL(indi,i,tt,gen) + 
     &                 (2*p1k2*(muu3**2*(-6 + x1) - (-2 + x1)*(1 + muu4**2 - x2)) - 2*p1k1*x2*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*((-1 + muu1**2 + muu2**2)*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*x1 - muu3**2*(6*muu4**2 + x1)) + 
     &                       (-2 + muu1**2*(6*muu3**2 - 2*muu4**2 - x1) - muu2**2*(-2 + x1) + x1)*x2 + (1 + muu1**2 - muu2**2)*x2**2))*GR(indi,i,tt,gen))*HR(indj,i,tt,gen)) + 
     &           2*dconjg(DL(j))*(muu4*GR(indi,i,tt,gen)*(-(muu2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*
     &                    HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HR(indj,i,tt,gen)) + 
     &              muu1*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*HL(indj,i,tt,gen) + 
     &                 muu3*(1 + muu3**2 + muu4**2 - x2)*(-2*p1k2 + sqrtS**2*x2)*HR(indj,i,tt,gen)))) + 
     &        dconjg(CL(j))*(dconjg(DL(j))*(GL(indi,i,tt,gen)*((2*p1k2*(muu3**2*(-6 + x1) - (-2 + x1)*(1 + muu4**2 - x2)) - 2*p1k1*x2*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*((-1 + muu1**2 + muu2**2)*(-1 + muu3**4 + muu4**4 + x1 + muu4**2*x1 - muu3**2*(6*muu4**2 + x1)) + 
     &                       (-2 + muu1**2*(6*muu3**2 - 2*muu4**2 - x1) - muu2**2*(-2 + x1) + x1)*x2 + (1 + muu1**2 - muu2**2)*x2**2))*HL(indj,i,tt,gen) - 
     &                 4*muu2*muu3*p1k2*(1 + muu3**2 + muu4**2 - x2)*HR(indj,i,tt,gen)) + 
     &              2*muu1*muu4*GR(indi,i,tt,gen)*(((1 - muu1**2 + muu2**2)*sqrtS**2 + 2*p1k2*(1 - 3*muu3**2 + muu4**2 - x2) - 2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + 
     &                    sqrtS**2*(muu4**2 - x1 - 2*x2 - (muu4**2 - x2)*(muu1**2 - muu2**2 + x1 + x2) + muu3**2*(-3 + 3*muu1**2 - 3*muu2**2 + x1 + 3*x2)))*HL(indj,i,tt,gen) - 
     &                 2*muu2*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HR(indj,i,tt,gen))) + 
     &           2*dconjg(DR(j))*(muu1*GR(indi,i,tt,gen)*(muu3*(1 + muu3**2 + muu4**2 - x2)*(-2*p1k2 + sqrtS**2*x2)*HL(indj,i,tt,gen) + 
     &                 muu2*sqrtS**2*(muu3**4 - muu3**2*(6*muu4**2 + x1 - 3*x2) + (-1 + muu4**2 + x1)*(1 + muu4**2 - x2))*HR(indj,i,tt,gen)) + 
     &              muu4*GL(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*(1 + muu3**2 + muu4**2 - x2)*HL(indj,i,tt,gen)) - 
     &                 muu2*(2*p1k1*(-1 + muu3**2 - muu4**2 + x2) + (2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(-1 + 3*muu3**2 - muu4**2 + x2))*HR(indj,i,tt,gen))))))/
     &    ((-1 + muu3**2 - muu4**2 + x2)**2*dreal(DenSmH(j))*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))))   )

** USchannel bremsstrahlung 14

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   (CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DR(j))*(-(muu4*GL(indi,i,tt,gen)*
     &               (2*muu3*(2*p1k1 + 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 - muu1**2*(x1 + x2)))*HL(indj,i,tt,gen) + 
     &                 muu2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 + muu3**2*(4*p1k1 + 12*p1k2 - 5*(1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                    muu4**2*(-4*(p1k1 + p1k2) + (1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 2*p1k1*x2 + 
     &                    (4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*x2)*HR(indj,i,tt,gen))) + 
     &            muu1*GR(indi,i,tt,gen)*(muu3*(2*p1k2*x1 - muu4**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x2)) - 2*p1k1*x2 + (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2) + 
     &                  muu3**2*(-4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*HR(indj,i,tt,gen))) + 
     &         dconjg(DL(j))*(muu1*muu4*GR(indi,i,tt,gen)*((4*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 + (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2)*x1 + 
     &                  muu4**2*(4*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) - (2*p1k1 + 4*p1k2 + (3 - muu1**2 + muu2**2)*sqrtS**2)*x2 + 2*sqrtS**2*x2**2 + 
     &                  muu3**2*(-4*p1k1 - 12*p1k2 + sqrtS**2*(-5 + 5*muu1**2 - 5*muu2**2 + 2*x1 + 6*x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + 4*p1k2*x1 - 2*muu1**2*sqrtS**2*x1*x2 + 
     &                  muu4**2*(-2*p1k2*x1 + 2*p1k1*(-2 + x2) + sqrtS**2*(-(muu1**2*(-1 + 6*muu3**2 - 3*x1 + x2)) - (-1 + muu2**2)*(-1 + 6*muu3**2 - x1 + x2))) + 
     &                  muu3**2*(2*p1k2*(-6 + x1) - 2*p1k1*x2 + sqrtS**2*(-((-1 + muu2**2)*(-1 + x1 + x2)) + muu1**2*(1 - x1 + 5*x2))))*HL(indj,i,tt,gen) + 
     &               muu2*muu3*(muu4**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu3**2*(-4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) + 2*p1k2*x1 - 2*p1k1*x2 + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HR(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(-(dconjg(DL(j))*(muu1*GL(indi,i,tt,gen)*(-2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*
     &                  HL(indj,i,tt,gen) + muu3*(-2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*x2)) + muu4**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x2)) + 
     &                    2*p1k1*x2 - (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HR(indj,i,tt,gen)) + 
     &              muu4*GR(indi,i,tt,gen)*(muu2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 + muu3**2*(4*p1k1 + 12*p1k2 - 5*(1 + muu1**2 - muu2**2)*sqrtS**2) + 
     &                    muu4**2*(-4*(p1k1 + p1k2) + (1 + muu1**2 - muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 2*p1k1*x2 + 
     &                    (4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*x2)*HL(indj,i,tt,gen) + 
     &                 2*muu3*(2*p1k1 + 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 - muu1**2*(x1 + x2)))*HR(indj,i,tt,gen)))) + 
     &         dconjg(DR(j))*(muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HL(indj,i,tt,gen) + 
     &               (4*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 + (2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2)*x1 + 
     &                  muu4**2*(4*p1k1 + 4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) - (2*p1k1 + 4*p1k2 + (3 - muu1**2 + muu2**2)*sqrtS**2)*x2 + 2*sqrtS**2*x2**2 + 
     &                  muu3**2*(-4*p1k1 - 12*p1k2 + sqrtS**2*(-5 + 5*muu1**2 - 5*muu2**2 + 2*x1 + 6*x2)))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(muu2*muu3*(muu4**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu3**2*(-4*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2) + 2*p1k2*x1 - 2*p1k1*x2 + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HL(indj,i,tt,gen) + 
     &               ((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + (-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + 4*p1k2*x1 - 2*muu1**2*sqrtS**2*x1*x2 + 
     &                  muu4**2*(-2*p1k2*x1 + 2*p1k1*(-2 + x2) + sqrtS**2*(-(muu1**2*(-1 + 6*muu3**2 - 3*x1 + x2)) - (-1 + muu2**2)*(-1 + 6*muu3**2 - x1 + x2))) + 
     &                  muu3**2*(2*p1k2*(-6 + x1) - 2*p1k1*x2 + sqrtS**2*(-((-1 + muu2**2)*(-1 + x1 + x2)) + muu1**2*(1 - x1 + 5*x2))))*HR(indj,i,tt,gen))))))/
     &  ((1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmH(j))*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))    )

** USchannel bremsstrahlung 23

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   (CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(muu1*muu4*GR(indi,i,tt,gen)*
     &             ((muu4**2*(4*p1k1 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*x1)) + muu3**2*(4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1)) + 2*p1k2*x1 - 2*p1k1*x2 - 
     &                  (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(((-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + muu3**2*((-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 - x1) - 2*p1k2*(-2 + x1)) + 
     &                  (2*(-2 + muu3**2)*p1k1 + sqrtS**2*((-3 + muu1**2 + 3*muu2**2)*muu3**2 - 2*(-1 + muu2**2)*x1))*x2 - 
     &                  muu4**2*(-2*p1k2*x1 + 2*p1k1*(-6 + x2) + sqrtS**2*((-1 + muu2**2)*(-1 + 6*muu3**2 - 5*x1 + x2) + muu1**2*(-1 + 6*muu3**2 + x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*muu3*(4*muu4**2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 - 5*muu4**2*sqrtS**2 - 5*muu1**2*muu4**2*sqrtS**2 + 5*muu2**2*muu4**2*sqrtS**2 + 
     &                  muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*x2 - 
     &                  2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2))*HR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(muu1*GR(indi,i,tt,gen)*(-(muu3*(-4*muu4**2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 - 5*muu4**2*sqrtS**2 + 5*muu1**2*muu4**2*sqrtS**2 - 
     &                    5*muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 3*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 6*muu4**2*sqrtS**2*x1 + 2*sqrtS**2*x1**2 + 
     &                    muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) + (1 - muu1**2 + muu2**2 + 2*muu4**2)*sqrtS**2*x2 + 2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2)
     &                    )*HL(indj,i,tt,gen)) + 2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*HR(indj,i,tt,gen)) + 
     &            muu4*GL(indi,i,tt,gen)*(-2*muu3*(-2*p1k1 - 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 + x1 + x2 - muu2**2*(x1 + x2)))*
     &                HL(indj,i,tt,gen) + muu2*(muu4**2*(4*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu3**2*(4*p1k1 + (-1 - muu1**2 + muu2**2)*sqrtS**2) + 2*p1k2*x1 - 2*p1k1*x2 - 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HR(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(dconjg(DL(j))*(muu1*GL(indi,i,tt,gen)*(2*muu2*sqrtS**2*(muu3**4 + muu4**4 - muu3**2*(-1 + 6*muu4**2 + x1 - 2*x2) + muu4**2*(1 + 2*x1 - x2) - x1*x2)*
     &                HL(indj,i,tt,gen) - muu3*(-4*muu4**2*p1k2 + (1 - muu1**2 + muu2**2)*sqrtS**2 - 5*muu4**2*sqrtS**2 + 5*muu1**2*muu4**2*sqrtS**2 - 5*muu2**2*muu4**2*sqrtS**2 - 
     &                  2*p1k2*x1 - 3*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 6*muu4**2*sqrtS**2*x1 + 2*sqrtS**2*x1**2 + 
     &                  muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1 - 2*x2)) + (1 - muu1**2 + muu2**2 + 2*muu4**2)*sqrtS**2*x2 + 2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2))*
     &                HR(indj,i,tt,gen)) + muu4*GR(indi,i,tt,gen)*(muu2*(muu4**2*(4*p1k1 + (1 + muu1**2 - muu2**2)*sqrtS**2) + muu3**2*(4*p1k1 + (-1 - muu1**2 + muu2**2)*sqrtS**2) + 
     &                  2*p1k2*x1 - 2*p1k1*x2 - (1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HL(indj,i,tt,gen) - 
     &               2*muu3*(-2*p1k1 - 2*p1k2 + sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**2 + (-1 + muu1**2 + muu2**2)*muu4**2 + x1 + x2 - muu2**2*(x1 + x2)))*HR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(-2*muu3**2 - 2*muu4**2 + x1 + x2)*HL(indj,i,tt,gen) + 
     &               (muu4**2*(4*p1k1 + sqrtS**2*(-1 + muu1**2 - muu2**2 - 2*x1)) + muu3**2*(4*p1k1 + sqrtS**2*(1 - muu1**2 + muu2**2 - 2*x1)) + 2*p1k2*x1 - 2*p1k1*x2 - 
     &                  (-1 + muu1**2 - muu2**2)*sqrtS**2*(-1 + x1 + x2))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(muu2*muu3*(4*muu4**2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2 - 5*muu4**2*sqrtS**2 - 5*muu1**2*muu4**2*sqrtS**2 + 5*muu2**2*muu4**2*sqrtS**2 + 
     &                  muu3**2*(-4*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + (1 + muu1**2 - muu2**2)*sqrtS**2*x2 - 
     &                  2*p1k1*(2 + 2*muu3**2 - 6*muu4**2 - 2*x1 + x2))*HL(indj,i,tt,gen) + 
     &               ((-1 + muu1**2 + muu2**2)*muu4**4*sqrtS**2 + muu3**2*((-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 - x1) - 2*p1k2*(-2 + x1)) + 
     &                  (2*(-2 + muu3**2)*p1k1 + sqrtS**2*((-3 + muu1**2 + 3*muu2**2)*muu3**2 - 2*(-1 + muu2**2)*x1))*x2 - 
     &                  muu4**2*(-2*p1k2*x1 + 2*p1k1*(-6 + x2) + sqrtS**2*((-1 + muu2**2)*(-1 + 6*muu3**2 - 5*x1 + x2) + muu1**2*(-1 + 6*muu3**2 + x1 + x2))))*HR(indj,i,tt,gen))))))/
     &  ((1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmH(j))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))  )

** USchannel bremsstrahlung 24

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   -((CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(-2*muu1*muu4*(1 + muu3**2 + muu4**2 - x1)*GR(indi,i,tt,gen)*
     &               ((-2*p1k1 + sqrtS**2*x1)*HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &              GL(indi,i,tt,gen)*((2*p1k2*x1*(-1 - muu3**2 + muu4**2 + x1) + 2*p1k1*(-2 + 6*muu4**2 + 2*x1 + muu3**2*(-2 + x2) + x2 - (muu4**2 + x1)*x2) + 
     &                    sqrtS**2*(muu1**2*(muu3**4 + (-1 + muu4**2 + x1)*(1 + muu4**2 - x1 - x2) + muu3**2*(-6*muu4**2 + x2)) + 
     &                       (-1 + muu2**2)*(-1 + muu3**4 + muu4**4 + 6*muu4**2*x1 + x1**2 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - 2*x1 + x2))))*HL(indj,i,tt,gen) + 
     &                 2*muu2*muu3*((1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(-1 - muu3**2 + muu4**2 + x1) + 2*p1k1*(-1 - muu3**2 + 3*muu4**2 + x1))*
     &                  HR(indj,i,tt,gen))) + 2*dconjg(DR(j))*(muu4*(1 + muu3**2 + muu4**2 - x1)*GL(indi,i,tt,gen)*
     &               (-((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HL(indj,i,tt,gen)) + 2*muu2*p1k1*HR(indj,i,tt,gen)) + 
     &              muu1*GR(indi,i,tt,gen)*(-(muu3*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                      sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*HL(indj,i,tt,gen)) + 
     &                 muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*HR(indj,i,tt,gen)))) + 
     &        dconjg(CR(j))*(2*dconjg(DL(j))*(muu4*(1 + muu3**2 + muu4**2 - x1)*GR(indi,i,tt,gen)*
     &               (2*muu2*p1k1*HL(indj,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &              muu1*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(-1 + muu3**4 + muu4**4 + x1 + 3*muu4**2*x1 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - x1 + x2))*HL(indj,i,tt,gen) - 
     &                 muu3*(2*p1k1*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(1 + muu3**2 - muu4**2 - x1) + 
     &                    sqrtS**2*(-((1 + muu3**2 - 3*muu4**2 - x1)*(-1 + muu1**2 - muu2**2 + x1)) + (-1 - muu3**2 + muu4**2 + x1)*x2))*HR(indj,i,tt,gen))) + 
     &           dconjg(DR(j))*(-2*muu1*muu4*(1 + muu3**2 + muu4**2 - x1)*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*HL(indj,i,tt,gen) + (-2*p1k1 + sqrtS**2*x1)*HR(indj,i,tt,gen)) + 
     &              GR(indi,i,tt,gen)*(2*muu2*muu3*((1 + muu1**2 - muu2**2)*sqrtS**2*(1 + muu3**2 - 3*muu4**2 - x1) + 2*p1k2*(-1 - muu3**2 + muu4**2 + x1) + 
     &                    2*p1k1*(-1 - muu3**2 + 3*muu4**2 + x1))*HL(indj,i,tt,gen) + 
     &                 (2*p1k2*x1*(-1 - muu3**2 + muu4**2 + x1) + 2*p1k1*(-2 + 6*muu4**2 + 2*x1 + muu3**2*(-2 + x2) + x2 - (muu4**2 + x1)*x2) + 
     &                    sqrtS**2*(muu1**2*(muu3**4 + (-1 + muu4**2 + x1)*(1 + muu4**2 - x1 - x2) + muu3**2*(-6*muu4**2 + x2)) + 
     &                       (-1 + muu2**2)*(-1 + muu3**4 + muu4**4 + 6*muu4**2*x1 + x1**2 + x2 - (muu4**2 + x1)*x2 + muu3**2*(-6*muu4**2 - 2*x1 + x2))))*HR(indj,i,tt,gen))))))/
     &    ((-1 - muu3**2 + muu4**2 + x1)**2*dreal(DenSmH(j))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))))  )

** USchannel bremsstrahlung 33

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(   -(CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(GL(indi,i,tt,gen)
     &                   *((2*p1k1*(8*p1k2 + sqrtS**2*(4*(-1 + muu1**2)*muu4**2 - (1 + 4*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*x2 + x2**2)) + 
     &                   sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-3 + 2*muu2**2 + muu4**2 - x2) + 
     &                      muu3**2*(2*p1k2*(4 - 4*muu2**2 + x1) + sqrtS**2*(muu1**2*(2 + 2*muu2**2 + 6*muu4**2 - x1 - 6*x2) - (-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + x1 + 2*x2))) + 
     &                      sqrtS**2*(muu1**2*(1 + muu4**4 - x1 + muu4**2*(2 - 5*x1 - 2*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 2*x2 + 5*x1*x2 + x2**2) - 
     &                         (-1 + muu2)*(1 + muu2)*(-muu4**4 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - (-1 + x2)*(-1 + x1 + x2) + muu4**2*(-2 + 5*x1 + 2*x2)))))*HL(indj,i,tt,gen) + 
     &                4*muu2*muu3*(p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 + muu3**2 + muu4**2 - x2)))*
     &                 HR(indj,i,tt,gen)) + 2*muu1*muu4*GR(indi,i,tt,gen)*((-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &                   sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*
     &                 HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*HR(indj,i,tt,gen))) + 
     &          2*dconjg(DR(j))*(muu4*GL(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*(2*muu2**2*p1k2 + sqrtS**2 - p1k2*x1) + 
     &                   muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) - 2*p1k1*x2 + 
     &                   sqrtS**2*(-2*muu2**4 - 2*muu3**2 - 2*muu4**2 + x1 + muu2**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2) + 2*x2))*HL(indj,i,tt,gen) + 
     &                muu2*(8*p1k1*p1k2 - 2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k1*sqrtS**2 + 
     &                   sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                 HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                   sqrtS**2*(2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + 
     &                      sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - (muu3**2 + muu4**2)*x2 + x2**2 + 
     &                         muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &                muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 + muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + 
     &                   muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x2)*(-1 + x1 + x2))*HR(indj,i,tt,gen)))) + 
     &       dconjg(CR(j))*(2*dconjg(DL(j))*(muu4*GR(indi,i,tt,gen)*(muu2*(8*p1k1*p1k2 - 2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k1*sqrtS**2 + 
     &                   sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                 HL(indj,i,tt,gen) + muu3*sqrtS**2*(-2*(2*muu2**2*p1k2 + sqrtS**2 - p1k2*x1) + muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) - 
     &                   2*p1k1*x2 + sqrtS**2*(-2*muu2**4 - 2*muu3**2 - 2*muu4**2 + x1 + muu2**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2) + 2*x2))*HR(indj,i,tt,gen)) + 
     &             muu1*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 + muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + 
     &                   muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - sqrtS**2*(muu1**2 - muu2**2 - x2)*(-1 + x1 + x2))*HL(indj,i,tt,gen) + 
     &                muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + sqrtS**2*
     &                    (2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + sqrtS**2*
     &                       (-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - (muu3**2 + muu4**2)*x2 + x2**2 + 
     &                         muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HR(indj,i,tt,gen))) + 
     &          dconjg(DR(j))*(2*muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*HL(indj,i,tt,gen) + 
     &                (-2*p1k1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &                   sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-1 + 2*muu1**2 + muu4**2 + x2))))*
     &                 HR(indj,i,tt,gen)) + GR(indi,i,tt,gen)*(4*muu2*muu3*(p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + 
     &                   sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 + muu3**2 + muu4**2 - x2)))*HL(indj,i,tt,gen) + 
     &                (2*p1k1*(8*p1k2 + sqrtS**2*(4*(-1 + muu1**2)*muu4**2 - (1 + 4*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*x2 + x2**2)) + 
     &                   sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-3 + 2*muu2**2 + muu4**2 - x2) + 
     &                      muu3**2*(2*p1k2*(4 - 4*muu2**2 + x1) + sqrtS**2*(muu1**2*(2 + 2*muu2**2 + 6*muu4**2 - x1 - 6*x2) - (-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + x1 + 2*x2))) + 
     &                      sqrtS**2*(muu1**2*(1 + muu4**4 - x1 + muu4**2*(2 - 5*x1 - 2*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 2*x2 + 5*x1*x2 + x2**2) - 
     &                         (-1 + muu2)*(1 + muu2)*(-muu4**4 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - (-1 + x2)*(-1 + x1 + x2) + muu4**2*(-2 + 5*x1 + 2*x2)))))*HR(indj,i,tt,gen))))))/
     &  (2.*(-1 + muu3**2 - muu4**2 + x2)*dreal(DenSmH(j))*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

** USchannel bremsstrahlung 34

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS + 2*dreal(  (CF*gs**2*NC*(dconjg(CL(j))*(dconjg(DL(j))*(GL(indi,i,tt,gen)
     &                  *((8*p1k1*(2*p1k2 + (-1 + muu1**2)*muu4**2*sqrtS**2) - 2*p1k1*sqrtS**2*(1 + 2*muu1**2 + muu3**2 + muu4**2 - x1)*x2 + 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 - 2*p1k2*x1*(3 + 2*muu1**2 - 4*muu2**2 - muu4**2 + x1) + 
     &                     muu3**2*(2*p1k2*(4 - 4*muu2**2 + x1) + sqrtS**2*(-2*muu1**4 + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - 2*x1 - 5*x2) + (-1 + muu2**2)*(2 + 6*muu4**2 - 2*x1 - x2))) + 
     &                     sqrtS**2*(2*muu1**4*(-1 - muu4**2 + x1 + x2) + (-1 + muu2)*(1 + muu2)*(muu4**4 + (-1 + x1)*(-1 + x1 + x2) - muu4**2*(-2 + 6*x1 + x2)) + 
     &                        muu1**2*(3 + muu4**4 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 3*x2 - muu4**2*(-4 + 6*x1 + x2) + x1*(-4 + x1 + 5*x2)))))*HL(indj,i,tt,gen) + 
     &               4*muu2*muu3*(p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + sqrtS**2*
     &                   (sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*HR(indj,i,tt,gen)) + 
     &            2*muu1*muu4*GR(indi,i,tt,gen)*((2*p1k1*(4*p1k2 - sqrtS**2*(-1 + muu3**2 + muu4**2 - x1 + 2*x2)) + 
     &                  sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x2) + x1*(-1 + muu4**2 - x1 + 2*x2))))*
     &                HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*HR(indj,i,tt,gen))) + 
     &         2*dconjg(DR(j))*(muu4*GL(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*muu1**4*sqrtS**2 + 2*p1k2*x1 + muu1**2*(4*p1k1 + sqrtS**2*(2*(2 + muu2**2 + muu3**2 + muu4**2 - 2*x1) - x2)) + 
     &                  muu2**2*(-4*p1k2 + sqrtS**2*(2*(1 + muu3**2 + muu4**2 - x1) - x2)) - 2*p1k1*x2 + sqrtS**2*(-2*(1 + muu3**2 + muu4**2 - x1) + x2))*HL(indj,i,tt,gen) + 
     &               muu2*(sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*
     &                HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                  sqrtS**2*(2*p1k2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1) + 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - (muu3**2 + muu4**2 - x1)*x2 + 
     &                        muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x1)*(-1 + x1 + x2) + 
     &                  muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HR(indj,i,tt,gen))))
     &        + dconjg(CR(j))*(-2*dconjg(DL(j))*(muu4*GR(indi,i,tt,gen)*(muu2*
     &                (sqrtS**2*(muu3**2*(-4*p1k2 + (1 + 3*muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) + 
     &                  2*p1k1*(-4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*HL(indj,i,tt,gen) + 
     &               muu3*sqrtS**2*(2*muu1**4*sqrtS**2 - 2*p1k2*x1 + sqrtS**2*(2*(1 + muu3**2 + muu4**2 - x1) - x2) + 2*p1k1*x2 + 
     &                  muu1**2*(-4*p1k1 + sqrtS**2*(-2*(2 + muu2**2 + muu3**2 + muu4**2 - 2*x1) + x2)) + muu2**2*(4*p1k2 + sqrtS**2*(-2*(1 + muu3**2 + muu4**2 - x1) + x2)))*
     &                HR(indj,i,tt,gen)) + muu1*GL(indi,i,tt,gen)*(-(muu2*sqrtS**2*
     &                  (muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - 2*p1k1*x2 + sqrtS**2*(muu1**2 - muu2**2 + x1)*(-1 + x1 + x2) + 
     &                    muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HL(indj,i,tt,gen))
     &                 + muu3*(p1k1*(-8*p1k2 + 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                  sqrtS**2*(p1k2*(4*muu1**2 - 2*(1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1)) - 
     &                     sqrtS**2*(-1 - muu3**2 + x1 + muu4**2*(-1 + 2*x1) + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - (muu3**2 + muu4**2 - x1)*x2 + 
     &                        muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(2*muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*HL(indj,i,tt,gen) + 
     &               (2*p1k1*(4*p1k2 - sqrtS**2*(-1 + muu3**2 + muu4**2 - x1 + 2*x2)) + 
     &                  sqrtS**2*(-4*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + sqrtS**2*(2*muu2**2*(1 + muu4**2 - x2) + x1*(-1 + muu4**2 - x1 + 2*x2))))*
     &                HR(indj,i,tt,gen)) + GR(indi,i,tt,gen)*(4*muu2*muu3*(p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + 
     &                  sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*HL(indj,i,tt,gen) + 
     &               (8*p1k1*(2*p1k2 + (-1 + muu1**2)*muu4**2*sqrtS**2) - 2*p1k1*sqrtS**2*(1 + 2*muu1**2 + muu3**2 + muu4**2 - x1)*x2 + 
     &                  sqrtS**2*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 - 2*p1k2*x1*(3 + 2*muu1**2 - 4*muu2**2 - muu4**2 + x1) + 
     &                     muu3**2*(2*p1k2*(4 - 4*muu2**2 + x1) + sqrtS**2*(-2*muu1**4 + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - 2*x1 - 5*x2) + (-1 + muu2**2)*(2 + 6*muu4**2 - 2*x1 - x2))) + 
     &                     sqrtS**2*(2*muu1**4*(-1 - muu4**2 + x1 + x2) + (-1 + muu2)*(1 + muu2)*(muu4**4 + (-1 + x1)*(-1 + x1 + x2) - muu4**2*(-2 + 6*x1 + x2)) + 
     &                        muu1**2*(3 + muu4**4 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - 3*x2 - muu4**2*(-4 + 6*x1 + x2) + x1*(-4 + x1 + 5*x2)))))*HR(indj,i,tt,gen))))))/
     &  (2.*(1 + muu3**2 - muu4**2 - x1)*dreal(DenSmH(j))*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

        if (cutoff.lt.zerocheck) then

** USchannel dipole with H 312

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal(  (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) - 
     &           4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &           (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &         (dconjg(CL(j))*dconjg(DL(j))*(((-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                   (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                      x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*GL(indi,i,tt,gen)
     &                 )/(4.*sqrtS**4*(4*muu4**2 - x2**2)) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*HL(indj,i,tt,gen) + 
     &           dconjg(CR(j))*dconjg(DR(j))*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + 
     &              ((-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                   (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                      x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*GR(indi,i,tt,gen)
     &                 )/(4.*sqrtS**4*(4*muu4**2 - x2**2)))*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu4**2 - x2**2)) - 
     &      ((4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &           (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &              x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*
     &         (muu2*muu4*sqrtS**2*(dconjg(CR(j))*dconjg(DL(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DR(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)) - 
     &           ((4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 - 
     &                (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) - 
     &                4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &              (dconjg(CL(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(4.*sqrtS**4*(4*muu4**2 - x2**2))))
     &        /(sqrtS**4*(4*muu4**2 - x2**2)) + sqrtS**2*((muu2*muu3*(-4*p1k2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                 x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*
     &            (dconjg(CR(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu4**2 - x2**2)) - 
     &         (muu1*muu3*(4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 + 
     &              (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) + 
     &              4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &              (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*
     &            (dconjg(CL(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu4**2 - x2**2)) + 
     &         sqrtS**2*(dconjg(CL(j))*((-1 + muu1**2 + muu2**2)*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &               2*muu1*muu2*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen)) + 
     &            dconjg(CR(j))*(2*muu1*muu2*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &               (-1 + muu1**2 + muu2**2)*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*dreal(DenSmH(j))*(-((muu2**2 + muu3**2)*sqrtS**2) + 
     &      (4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 - (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) - 
     &         4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))/
     &       (2.*sqrtS**4*(4*muu4**2 - x2**2)) + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )
 
** USchannel dipole with H 321

	ChiChi2QQ_RadiationUS = ChiChi2QQ_RadiationUS - 2*dreal(  (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) + (Sqrt((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*
     &         (-3*muu3**4 + muu4**4 + muu4**2*(x1 - x2) - (-1 + x1)*(-1 + 2*x1 + x2) + muu3**2*(-4 - 2*muu4**2 + 5*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (((4*(1 - muu1**2 + muu2**2)*muu3**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 + 
     &           4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &           (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*
     &         (dconjg(CL(j))*dconjg(DL(j))*(((4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                   (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                      x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*GL(indi,i,tt,gen))
     &                /(4.*sqrtS**4*(4*muu3**2 - x1**2)) + muu1*muu4*sqrtS**2*GR(indi,i,tt,gen))*HL(indj,i,tt,gen) + 
     &           dconjg(CR(j))*dconjg(DR(j))*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + 
     &              ((4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                   (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                      x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*GR(indi,i,tt,gen))
     &                /(4.*sqrtS**4*(4*muu3**2 - x1**2)))*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu3**2 - x1**2)) - 
     &      4*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))*
     &       (muu2*muu4*sqrtS**2*(dconjg(CR(j))*dconjg(DL(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DR(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)) + 
     &         ((4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &              4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &              (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*
     &            (dconjg(CL(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(4.*sqrtS**4*(4*muu3**2 - x1**2))) + 
     &      sqrtS**2*((muu2*muu3*(4*p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                 x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*
     &            (dconjg(CR(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CL(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu3**2 - x1**2)) - 
     &         (muu1*muu3*(4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 - 
     &              4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &              (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &              (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*
     &            (dconjg(CL(j))*dconjg(DR(j))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + dconjg(CR(j))*dconjg(DL(j))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen)))/(sqrtS**4*(4*muu3**2 - x1**2)) + 
     &         sqrtS**2*(dconjg(CL(j))*((-1 + muu1**2 + muu2**2)*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &               2*muu1*muu2*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen)) + 
     &            dconjg(CR(j))*(2*muu1*muu2*((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &               (-1 + muu1**2 + muu2**2)*(2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen))))))/
     &  (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*dreal(DenSmH(j))*(-((muu2**2 + muu3**2)*sqrtS**2) + 
     &      (4*(1 - muu1**2 + muu2**2)*muu3**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 + 
     &         4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &         (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &         (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen))))   )

          endif


	  enddo
	enddo

	end

************************************************

	double precision function ChiChi2QQ_RadiationTU()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j

	double precision temp1,temp2
	

	ChiChi2QQ_RadiationTU = 0d0


***** itypeT,U follows the type of the sfermion exchanged in t,u-channel interactions
***** att gives the type of the antiquark final state
***** antitt is the opposite type to the final state quark type

	if (icha.eq.0) then
		itypeT = tt
		itypeU = tt
	else if (icha.eq.1) then
* tt is defined as normally as down-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt
		else if (tt.eq.4) then 
		      itypeT = tt
		      itypeU = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      itypeT = tt+1
		      itypeU = tt+1
		else if (tt.eq.4) then 
		      itypeT = tt-1
		      itypeU = tt-1
		endif	
	endif	

	
	do i = 1, 2
  	  do j = 1, 2

** TUchannel bremsstrahlung 11

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(  -((CF*gs**2*NC*(2*(-1 + muu1**2 + muu2**2)*muu3**2*sqrtS**4
     &        *(1 + muu3**2 - muu4**2 - x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu3**2*p1k2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu3**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        (-1 + muu1**2 + muu2**2)*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 4*muu3**2*p1k1*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu3**2*(-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + (-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*
     &         dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu3**2*p1k2*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*p1k2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 4*muu1*muu3**3*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*muu1*muu3*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu2*muu3**2*muu4*p1k1*sqrtS**2*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu2*muu3**2*muu4*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu2*muu4*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(1 - muu3**2 + muu4**2 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 8*muu1*muu2*muu3**3*muu4*sqrtS**4*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu3*muu4*sqrtS**4*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu1*muu3**2*muu4*sqrtS**2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu1*muu3**2*muu4*sqrtS**2*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 2*muu1*muu4*sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*
     &         dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*(-1 + muu1**2 + muu2**2)*muu3**3*muu4*sqrtS**4*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu1*muu2*muu3**2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu3**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu1*muu2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu2*muu3**3*p1k2*sqrtS**2*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu2*muu3*p1k2*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu2*muu3**3*p1k2*sqrtS**2*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu2*muu3*p1k2*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu1*muu2*muu3**2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu3**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu1*muu2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*(-1 + muu1**2 + muu2**2)*muu3**3*muu4*sqrtS**4*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu1*muu3**2*muu4*sqrtS**2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu1*muu3**2*muu4*sqrtS**2*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 2*muu1*muu4*sqrtS**2*(1 - muu3**2 + muu4**2 - x2)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*
     &         dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        8*muu1*muu2*muu3**3*muu4*sqrtS**4*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu3*muu4*sqrtS**4*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        8*muu2*muu3**2*muu4*p1k1*sqrtS**2*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu2*muu3**2*muu4*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu2*muu4*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(1 - muu3**2 + muu4**2 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 4*muu1*muu3**3*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*muu1*muu3*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*(-1 + muu1**2 + muu2**2)*muu3**2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu3**2*p1k2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu3**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        (-1 + muu1**2 + muu2**2)*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 4*muu3**2*p1k1*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu3**2*(-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + (-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*
     &         dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu3**2*p1k2*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*p1k2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen)))/(sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)**2*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &      (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))))   )

** TUchannel bremsstrahlung 12

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(   -((CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     &            (muu1*muu4*GR(indi,i,tt,gen)*(sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                    sqrtS**2*(3 - 3*muu3**2 + 3*muu4**2 - 3*x1 + 4*muu3**2*x1 - 5*x2 + 2*muu3**2*x2 - 2*muu4**2*x2 + 2*x1*x2 + 2*x2**2 + 
     &                       muu1**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2) - muu2**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2)))*HL(indj,i,tt,gen) + 
     &                 2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HR(indj,i,tt,gen)) + 
     &              GL(indi,i,tt,gen)*((-2*p1k1*(4*(-1 + muu3**2 + muu4**2)*p1k2 + sqrtS**2*((muu1**2 - muu2**2)*x2 - 2*muu4**2*(-2 + muu1**2 - muu2**2 + x2))) + 
     &                    sqrtS**2*(muu2**2*(2*p1k2*x1 - (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2)) + 
     &                       muu1**2*(-2*p1k2*x1 + sqrtS**2*(1 + muu4**2 - x1 + 2*muu4**2*x1 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2)) + muu1**4*sqrtS**2*(-1 - muu4**2 + x1 + x2) + 
     &                       muu3**2*(4*p1k2*(-2 + muu1**2 - muu2**2 + x1) + sqrtS**2*(-muu1**4 - (-1 + muu2**2)*(muu2**2 + 4*muu4**2) + muu1**2*(1 + 2*muu2**2 - 4*muu4**2 + 2*x2)))))*
     &                  HL(indj,i,tt,gen) + muu2*muu3*sqrtS**2*(sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                    muu4**2*(-8*p1k1 + 3*(1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*(-2 + x1) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - 
     &                    muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HR(indj,i,tt,gen))) + 
     &           dconjg(FR(indj,j,tt,gen))*(muu4*GL(indi,i,tt,gen)*(muu2*sqrtS**2*
     &                  (sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu4**2*sqrtS**2 + muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + 
     &                    muu3**2*(-8*p1k2 + 3*(1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 2*p1k1*(-2 + x2) - sqrtS**2*x2 - 
     &                    muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HL(indj,i,tt,gen) - 
     &                 2*muu3*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*(muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                    2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*HR(indj,i,tt,gen)) + 
     &              muu1*GR(indi,i,tt,gen)*(2*muu2*sqrtS**4*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))*HL(indj,i,tt,gen) + 
     &                 muu3*sqrtS**2*(-2*p1k2*(-2 + 4*muu4**2 + x1) + p1k1*(4*muu3**2 - 2*(-2 + 2*muu4**2 + 2*x1 + x2)) + 
     &                    sqrtS**2*(3 + 3*muu3**2 - 3*muu4**2 - 5*x1 - 2*muu3**2*x1 + 2*muu4**2*x1 + 2*x1**2 + muu2**2*(3 + 3*muu3**2 - 3*muu4**2 - 3*x1 - x2) - 3*x2 + 4*muu4**2*x2 + 
     &                       2*x1*x2 + muu1**2*(-3 - 3*muu3**2 + 3*muu4**2 + 3*x1 + x2)))*HR(indj,i,tt,gen)))) + 
     &        dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*(-(muu4*GR(indi,i,tt,gen)*
     &                 (2*muu3*(4*p1k1**2 + 4*p1k2**2 + sqrtS**4*(muu1**4 + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2) + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - x1 - x2)) - 
     &                      2*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) + 2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)))*HL(indj,i,tt,gen) - 
     &                   muu2*sqrtS**2*(sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu4**2*sqrtS**2 + muu1**2*muu4**2*sqrtS**2 - muu2**2*muu4**2*sqrtS**2 + 
     &                      muu3**2*(-8*p1k2 + 3*(1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 2*p1k1*(-2 + x2) - sqrtS**2*x2 - 
     &                      muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HR(indj,i,tt,gen))) + 
     &              muu1*GL(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*p1k2*(-2 + 4*muu4**2 + x1) + p1k1*(4*muu3**2 - 2*(-2 + 2*muu4**2 + 2*x1 + x2)) + 
     &                    sqrtS**2*(3 + 3*muu3**2 - 3*muu4**2 - 5*x1 - 2*muu3**2*x1 + 2*muu4**2*x1 + 2*x1**2 + muu2**2*(3 + 3*muu3**2 - 3*muu4**2 - 3*x1 - x2) - 3*x2 + 4*muu4**2*x2 + 
     &                       2*x1*x2 + muu1**2*(-3 - 3*muu3**2 + 3*muu4**2 + 3*x1 + x2)))*HL(indj,i,tt,gen) + 
     &                 2*muu2*sqrtS**4*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))*HR(indj,i,tt,gen))) + 
     &           dconjg(FR(indj,j,tt,gen))*(muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HL(indj,i,tt,gen) + 
     &                 sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                    sqrtS**2*(3 - 3*muu3**2 + 3*muu4**2 - 3*x1 + 4*muu3**2*x1 - 5*x2 + 2*muu3**2*x2 - 2*muu4**2*x2 + 2*x1*x2 + 2*x2**2 + 
     &                       muu1**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2) - muu2**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2)))*HR(indj,i,tt,gen)) + 
     &              GR(indi,i,tt,gen)*(muu2*muu3*sqrtS**2*(sqrtS**2 + muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**2*sqrtS**2 + muu1**2*muu3**2*sqrtS**2 - muu2**2*muu3**2*sqrtS**2 + 
     &                    muu4**2*(-8*p1k1 + 3*(1 + muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*(-2 + x1) - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 - 
     &                    muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2)*HL(indj,i,tt,gen) + 
     &                 (-2*p1k1*(4*(-1 + muu3**2 + muu4**2)*p1k2 + sqrtS**2*((muu1**2 - muu2**2)*x2 - 2*muu4**2*(-2 + muu1**2 - muu2**2 + x2))) + 
     &                    sqrtS**2*(muu2**2*(2*p1k2*x1 - (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2)) + 
     &                       muu1**2*(-2*p1k2*x1 + sqrtS**2*(1 + muu4**2 - x1 + 2*muu4**2*x1 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2)) + muu1**4*sqrtS**2*(-1 - muu4**2 + x1 + x2) + 
     &                       muu3**2*(4*p1k2*(-2 + muu1**2 - muu2**2 + x1) + sqrtS**2*(-muu1**4 - (-1 + muu2**2)*(muu2**2 + 4*muu4**2) + muu1**2*(1 + 2*muu2**2 - 4*muu4**2 + 2*x2)))))*
     &                  HR(indj,i,tt,gen))))))/
     &    (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &      (-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))))   )

** TUchannel bremsstrahlung 13

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(  (CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(2*sqrtS**2*dconjg(FR(indj,j,tt,gen))*
     &          (muu4*GL(indi,i,tt,gen)*(-2*muu2*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + 
     &                  muu1**2*sqrtS**4*(1 + muu4**2 - x2))*HL(indj,i,tt,gen) + 
     &               muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu1**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 - 2*p1k2*x1 + sqrtS**2*x1 + 
     &                  muu1**2*(4*p1k2 + sqrtS**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 4*x2)) + muu2**2*(-4*p1k1 + sqrtS**2*(2 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) + 2*p1k1*x2 + 
     &                  2*sqrtS**2*x2)*HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*
     &             (muu2*sqrtS**2*(-(muu1**2*sqrtS**2) + muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - 
     &                  sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + sqrtS**2*x2**2 + 
     &                  muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*HL(indj,i,tt,gen) + 
     &               muu3*(8*p1k2**2 + 2*p1k2*sqrtS**2*(3 - 2*muu1**2 + 2*muu2**2 + 3*muu3**2 + 3*muu4**2 - 2*x1 - 5*x2) + 
     &                  muu4**2*(4*p1k1*sqrtS**2 - sqrtS**4*(-2 + 2*muu1**2 + 2*x1 + 3*x2)) + sqrtS**4*(2*muu2**2*(1 + muu3**2 - x1 - x2) + x2*(-3 + 2*muu1**2 - 3*muu3**2 + 2*x1 + 3*x2)))*
     &                HR(indj,i,tt,gen))) + dconjg(FL(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*
     &             ((2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-1 - muu4**2 - muu4**2*x1 + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*HL(indj,i,tt,gen) - 
     &               2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**6 + 
     &                  2*p1k1*(16*p1k2**2 + 4*p1k2*sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 3*x2) - 
     &                     sqrtS**4*((1 - 2*muu1**2 + muu3**2 - x2)*x2 + muu4**2*(-4 + 4*muu2**2 + x2))) - 
     &                  muu3**2*(16*p1k2**2*sqrtS**2 + 2*p1k2*sqrtS**4*(-4*muu1**2 + 4*muu2**2 + x1 - 4*x2) + 
     &                     sqrtS**6*(2*muu1**4 - (-1 + muu2**2)*(2 + 6*muu4**2 - x1 - 2*x2) + muu1**2*(-4 - 2*muu2**2 - 6*muu4**2 + x1 + 6*x2))) + 
     &                  sqrtS**4*((-2*p1k2*x1 + (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*(1 + muu4**2 - x2) - 2*muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu1**2*(p1k2*(8 + 8*muu4**2 - 4*x1 - 8*x2) + sqrtS**2*(muu4**4 + muu4**2*(4 - 5*x1 - 6*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) + (-1 + x1 + x2)*(-3 + 5*x2)))))*
     &                HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) + sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                HR(indj,i,tt,gen)))) + dconjg(ER(indi,j,tt,gen))*(-2*sqrtS**2*dconjg(FL(indj,j,tt,gen))*
     &          (muu4*GR(indi,i,tt,gen)*(muu3*sqrtS**2*(2*sqrtS**2 + 2*muu1**4*sqrtS**2 + 2*muu3**2*sqrtS**2 + 2*muu4**2*sqrtS**2 + 2*p1k2*x1 - sqrtS**2*x1 - 2*p1k1*x2 - 2*sqrtS**2*x2 + 
     &                  muu2**2*(4*p1k1 + sqrtS**2*(-2 - 2*muu1**2 - 2*muu3**2 - 2*muu4**2 + x1 + 2*x2)) + muu1**2*(-4*p1k2 + sqrtS**2*(-4 - 2*muu3**2 - 2*muu4**2 + x1 + 4*x2)))*
     &                HL(indj,i,tt,gen) + 2*muu2*(muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + p1k1*(4*p1k2 + sqrtS**2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x2)) + 
     &                  muu1**2*sqrtS**4*(1 + muu4**2 - x2))*HR(indj,i,tt,gen)) + 
     &            muu1*GL(indi,i,tt,gen)*(muu3*(-8*p1k2**2 + 2*p1k2*sqrtS**2*(-3 + 2*muu1**2 - 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + 2*x1 + 5*x2) + 
     &                  sqrtS**4*(-2*muu2**2*(1 + muu3**2 - x1 - x2) + (3 - 2*muu1**2 + 3*muu3**2 - 2*x1 - 3*x2)*x2) + muu4**2*(-4*p1k1*sqrtS**2 + sqrtS**4*(-2 + 2*muu1**2 + 2*x1 + 3*x2)))*
     &                HL(indj,i,tt,gen) - muu2*sqrtS**2*(-(muu1**2*sqrtS**2) + muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 + muu1**2*sqrtS**2*x1 - 
     &                  muu2**2*sqrtS**2*x1 + 2*p1k1*x2 - sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + sqrtS**2*x2**2 + 
     &                  muu4**2*(-4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + x1 + 2*x2)) + muu3**2*(4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + 4*x2)))*HR(indj,i,tt,gen)))
     &          + dconjg(FR(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             (-2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x2))*HL(indj,i,tt,gen) + 
     &               (2*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x2)) + 
     &                  sqrtS**2*(-2*p1k2*x1 - muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) + 
     &                     sqrtS**2*(-1 - muu4**2 - muu4**2*x1 + muu2**2*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2) + x2 + x1*x2)))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*p1k2*(1 + muu3**2 + muu4**2 + x1 - x2) + sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*
     &                HL(indj,i,tt,gen) + ((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**6 + 
     &                  2*p1k1*(16*p1k2**2 + 4*p1k2*sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 3*x2) - 
     &                     sqrtS**4*((1 - 2*muu1**2 + muu3**2 - x2)*x2 + muu4**2*(-4 + 4*muu2**2 + x2))) - 
     &                  muu3**2*(16*p1k2**2*sqrtS**2 + 2*p1k2*sqrtS**4*(-4*muu1**2 + 4*muu2**2 + x1 - 4*x2) + 
     &                     sqrtS**6*(2*muu1**4 - (-1 + muu2**2)*(2 + 6*muu4**2 - x1 - 2*x2) + muu1**2*(-4 - 2*muu2**2 - 6*muu4**2 + x1 + 6*x2))) + 
     &                  sqrtS**4*((-2*p1k2*x1 + (-1 + muu2**2)*sqrtS**2*(1 + muu4**2 - x1 - x2))*(1 + muu4**2 - x2) - 2*muu1**4*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     muu1**2*(p1k2*(8 + 8*muu4**2 - 4*x1 - 8*x2) + sqrtS**2*(muu4**4 + muu4**2*(4 - 5*x1 - 6*x2) + 2*muu2**2*(1 + muu4**2 - x1 - x2) + (-1 + x1 + x2)*(-3 + 5*x2)))))*
     &                HR(indj,i,tt,gen))))))/
     &  (2.*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))   )

** TUchannel bremsstrahlung 21

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(    -((CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     &            (muu1*muu4*GR(indi,i,tt,gen)*(sqrtS**2*(sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu4**2*sqrtS**2 - muu1**2*muu4**2*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - 
     &                    2*p1k2*x1 - 3*sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + muu3**2*(8*p1k2 + sqrtS**2*(3 - 3*muu1**2 + 3*muu2**2 - 4*x2)) - 2*p1k1*(-2 + x2) - 
     &                    sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + 2*sqrtS**2*x1*x2)*HL(indj,i,tt,gen) + 
     &                 2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HR(indj,i,tt,gen)) + 
     &              GL(indi,i,tt,gen)*((p1k1*(-8*(-1 + muu3**2 + muu4**2)*p1k2 + 2*sqrtS**2*(2*(2 + muu1**2 - muu2**2)*muu4**2 + (-2 - muu1**2 + muu2**2 + 2*muu3**2)*x2)) + 
     &                    sqrtS**2*(muu2**2*sqrtS**2 - muu2**4*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - muu2**4*muu4**2*sqrtS**2 - 4*p1k2*x1 + 2*muu2**2*p1k2*x1 + 4*muu4**2*p1k2*x1 - 
     &                       muu2**2*sqrtS**2*x1 + muu2**4*sqrtS**2*x1 - 4*muu4**2*sqrtS**2*x1 + 4*muu2**2*muu4**2*sqrtS**2*x1 + 
     &                       muu3**2*(4*(2 + muu1**2 - muu2**2)*p1k2 + sqrtS**2*(-muu1**4 - (-1 + muu2**2)*(muu2**2 + 4*muu4**2 - 4*x2) + muu1**2*(1 + 2*muu2**2 - 4*muu4**2 - 2*x2))) - 
     &                       muu2**2*sqrtS**2*x2 + muu2**4*sqrtS**2*x2 + 2*sqrtS**2*x1*x2 - 2*muu2**2*sqrtS**2*x1*x2 + muu1**4*sqrtS**2*(-1 - muu4**2 + x1 + x2) + 
     &                       muu1**2*(-2*p1k2*x1 + sqrtS**2*(1 + muu4**2 - x1 - 2*muu4**2*x1 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2 + 2*x1*x2))))*HL(indj,i,tt,gen) + 
     &                 muu2*muu3*sqrtS**2*(2*p1k2*(-2 + 4*muu4**2 + x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(3 + 3*muu3**2 - 3*muu4**2 - 3*x1 - x2) + 
     &                    p1k1*(-4 - 4*muu3**2 + 4*muu4**2 + 4*x1 + 2*x2))*HR(indj,i,tt,gen))) + 
     &           dconjg(FR(indj,j,tt,gen))*(muu1*GR(indi,i,tt,gen)*(2*muu2*sqrtS**4*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))*HL(indj,i,tt,gen) + 
     &                 muu3*sqrtS**2*(sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu3**2*sqrtS**2 - muu1**2*muu3**2*sqrtS**2 + muu2**2*muu3**2*sqrtS**2 + 
     &                    muu4**2*(8*p1k1 + sqrtS**2*(3 - 3*muu1**2 + 3*muu2**2 - 4*x1)) - 2*p1k2*(-2 + x1) - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 - 2*p1k1*x2 - 
     &                    3*sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + 2*sqrtS**2*x1*x2)*HR(indj,i,tt,gen)) - 
     &              muu4*GL(indi,i,tt,gen)*(muu2*sqrtS**2*(-2*p1k1*(-2 + 4*muu3**2 + x2) - 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) + 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2))*HL(indj,i,tt,gen) + 
     &                 2*muu3*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                    sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*HR(indj,i,tt,gen)))) + 
     &        dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*(muu1*GL(indi,i,tt,gen)*
     &               (muu3*sqrtS**2*(sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu3**2*sqrtS**2 - muu1**2*muu3**2*sqrtS**2 + muu2**2*muu3**2*sqrtS**2 + 
     &                    muu4**2*(8*p1k1 + sqrtS**2*(3 - 3*muu1**2 + 3*muu2**2 - 4*x1)) - 2*p1k2*(-2 + x1) - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 - 2*p1k1*x2 - 
     &                    3*sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + 2*sqrtS**2*x1*x2)*HL(indj,i,tt,gen) + 
     &                 2*muu2*sqrtS**4*(1 - x1 + muu4**2*(1 + x1) - x2 + muu3**2*(1 - 4*muu4**2 + x2))*HR(indj,i,tt,gen)) + 
     &              muu4*GR(indi,i,tt,gen)*(-2*muu3*(4*p1k1**2 + 4*p1k2**2 - 2*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2) - 
     &                    2*p1k1*(-4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + x1 + x2)) + 
     &                    sqrtS**4*(muu1**4 + (-1 + muu2**2)*(1 + muu2**2 + muu3**2 + muu4**2 - 2*x1 - 2*x2) + muu1**2*(-2*muu2**2 + muu3**2 + muu4**2 + x1 + x2)))*HL(indj,i,tt,gen) + 
     &                 muu2*sqrtS**2*(2*p1k1*(-2 + 4*muu3**2 + x2) + 2*p1k2*(-2 + 2*muu3**2 - 2*muu4**2 + x1 + 2*x2) - 
     &                    (1 + muu1**2 - muu2**2)*sqrtS**2*(-3 + 3*muu3**2 - 3*muu4**2 + x1 + 3*x2))*HR(indj,i,tt,gen))) + 
     &           dconjg(FR(indj,j,tt,gen))*(muu1*muu4*GL(indi,i,tt,gen)*(2*muu2*muu3*sqrtS**4*(2*muu3**2 + 2*muu4**2 - x1 - x2)*HL(indj,i,tt,gen) + 
     &                 sqrtS**2*(sqrtS**2 - muu1**2*sqrtS**2 + muu2**2*sqrtS**2 + muu4**2*sqrtS**2 - muu1**2*muu4**2*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 3*sqrtS**2*x1 + 
     &                    muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + muu3**2*(8*p1k2 + sqrtS**2*(3 - 3*muu1**2 + 3*muu2**2 - 4*x2)) - 2*p1k1*(-2 + x2) - sqrtS**2*x2 + muu1**2*sqrtS**2*x2 - 
     &                    muu2**2*sqrtS**2*x2 + 2*sqrtS**2*x1*x2)*HR(indj,i,tt,gen)) + 
     &              GR(indi,i,tt,gen)*(muu2*muu3*sqrtS**2*(2*p1k2*(-2 + 4*muu4**2 + x1) + (1 + muu1**2 - muu2**2)*sqrtS**2*(3 + 3*muu3**2 - 3*muu4**2 - 3*x1 - x2) + 
     &                    p1k1*(-4 - 4*muu3**2 + 4*muu4**2 + 4*x1 + 2*x2))*HL(indj,i,tt,gen) + 
     &                 (p1k1*(-8*(-1 + muu3**2 + muu4**2)*p1k2 + 2*sqrtS**2*(2*(2 + muu1**2 - muu2**2)*muu4**2 + (-2 - muu1**2 + muu2**2 + 2*muu3**2)*x2)) + 
     &                    sqrtS**2*(muu2**2*sqrtS**2 - muu2**4*sqrtS**2 + muu2**2*muu4**2*sqrtS**2 - muu2**4*muu4**2*sqrtS**2 - 4*p1k2*x1 + 2*muu2**2*p1k2*x1 + 4*muu4**2*p1k2*x1 - 
     &                       muu2**2*sqrtS**2*x1 + muu2**4*sqrtS**2*x1 - 4*muu4**2*sqrtS**2*x1 + 4*muu2**2*muu4**2*sqrtS**2*x1 + 
     &                       muu3**2*(4*(2 + muu1**2 - muu2**2)*p1k2 + sqrtS**2*(-muu1**4 - (-1 + muu2**2)*(muu2**2 + 4*muu4**2 - 4*x2) + muu1**2*(1 + 2*muu2**2 - 4*muu4**2 - 2*x2))) - 
     &                       muu2**2*sqrtS**2*x2 + muu2**4*sqrtS**2*x2 + 2*sqrtS**2*x1*x2 - 2*muu2**2*sqrtS**2*x1*x2 + muu1**4*sqrtS**2*(-1 - muu4**2 + x1 + x2) + 
     &                       muu1**2*(-2*p1k2*x1 + sqrtS**2*(1 + muu4**2 - x1 - 2*muu4**2*x1 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2 + 2*x1*x2))))*HR(indj,i,tt,gen))))))/
     &    (sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &      (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))))   )

** TUchannel bremsstrahlung 22

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(   -((CF*gs**2*NC*(4*muu4**2*p1k2*(-2*p1k1 + sqrtS**2*x1)
     &        *dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu4**2*(-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) - (-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*
     &         dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*(-1 + muu1**2 + muu2**2)*muu4**2*sqrtS**4*(1 - muu3**2 + muu4**2 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu4**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        (-1 + muu1**2 + muu2**2)*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 4*muu4**2*p1k1*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu4**2*p1k1*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*p1k1*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 4*muu1*muu3*muu4**2*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu1*muu3*muu4**2*sqrtS**2*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 2*muu1*muu3*sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*
     &         dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu2*muu4**3*p1k1*sqrtS**2*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu2*muu4*p1k1*sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu1*muu2*muu3*muu4**3*sqrtS**4*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu1*muu2*muu3*muu4*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu1*muu4**3*sqrtS**2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*muu1*muu4*sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*(-1 + muu1**2 + muu2**2)*muu3*muu4**3*sqrtS**4*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        4*muu1*muu2*muu4**2*sqrtS**4*(1 - muu3**2 + muu4**2 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu4**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu1*muu2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        8*muu2*muu3*muu4**2*p1k2*sqrtS**2*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &        4*muu2*muu3*muu4**2*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &        2*muu2*muu3*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(1 + muu3**2 - muu4**2 - x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HL(indj,i,tt,gen) + 8*muu2*muu3*muu4**2*p1k2*sqrtS**2*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu2*muu3*muu4**2*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu2*muu3*sqrtS**2*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2)*sqrtS**2)*(1 + muu3**2 - muu4**2 - x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 4*muu1*muu2*muu4**2*sqrtS**4*(1 - muu3**2 + muu4**2 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu1*muu2*muu4**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu1*muu2*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*(-1 + muu1**2 + muu2**2)*muu3*muu4**3*sqrtS**4*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu1*muu4**3*sqrtS**2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*muu1*muu4*sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        8*muu1*muu2*muu3*muu4**3*sqrtS**4*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu1*muu2*muu3*muu4*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*dconjg(ELc(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        8*muu2*muu4**3*p1k1*sqrtS**2*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu2*muu4*p1k1*sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FL(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu1*muu3*muu4**2*sqrtS**2*(-2*p1k2 + sqrtS**2*x2)*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        4*muu1*muu3*muu4**2*sqrtS**2*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 2*muu1*muu3*sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ELc(indi,j,tt,gen))*
     &         dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu4**2*p1k2*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*muu4**2*(-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) - (-2*p1k1 - 2*p1k2 + (1 + muu1**2 - muu2**2)*sqrtS**2)*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 + sqrtS**2*x1)*dconjg(ER(indi,j,tt,gen))*
     &         dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*(-1 + muu1**2 + muu2**2)*muu4**2*sqrtS**4*(1 - muu3**2 + muu4**2 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &        2*(-1 + muu1**2 + muu2**2)*muu4**2*sqrtS**4*(1 + muu3**2 + muu4**2 - x1 - x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        (-1 + muu1**2 + muu2**2)*sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(-1 + muu3**2 - muu4**2 + x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen) + 4*muu4**2*p1k1*(-2*p1k2 + sqrtS**2*x2)*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        4*muu4**2*p1k1*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &        2*p1k1*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k1 - 2*p1k2 + sqrtS**2*(-1 + muu1**2 - muu2**2 + x1 + x2))*dconjg(ER(indi,j,tt,gen))*dconjg(FR(indj,j,tt,gen))*GR(indi,i,tt,gen)*
     &         HR(indj,i,tt,gen)))/(sqrtS**2*(-1 - muu3**2 + muu4**2 + x1)**2*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &      (2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))))  )

** TUchannel bremsstrahlung 23

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(   (CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(-2*sqrtS**2*dconjg(FR(indj,j,tt,gen))*
     &          (muu4*GL(indi,i,tt,gen)*(2*muu2*(4*p1k1**2 + muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + muu1**2*sqrtS**4*(1 + muu4**2 - x2) + 
     &                  p1k1*sqrtS**2*(-3 - 2*muu1**2 + 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + x1 + 2*x2))*HL(indj,i,tt,gen) + 
     &               muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu2**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 - 2*p1k2*x1 + 2*sqrtS**2*x1 + 
     &                  muu1**2*(4*p1k2 + sqrtS**2*(2 + 2*muu3**2 + 2*muu4**2 - 2*x1 - 3*x2)) + muu2**2*(-4*p1k1 + sqrtS**2*(4 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - 2*x1 - x2)) + 
     &                  2*p1k1*x2 + sqrtS**2*x2)*HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*
     &             (muu2*sqrtS**2*(muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 
     &                  sqrtS**2*x1**2 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 
     &                  2*p1k1*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2)*HL(indj,i,tt,gen) + 
     &               muu3*(-4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*
     &                   (2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) - 
     &                     sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + x1)*x2 + muu4**2*(2 - 2*muu1**2 - 2*x1 + x2))))*HR(indj,i,tt,gen))) + 
     &         dconjg(FL(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*
     &             ((-2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                  sqrtS**2*(2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) - 
     &                     sqrtS**2*(-1 - muu4**2 - muu4**2*x1 + x1**2 + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((16*p1k1**2*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) - 
     &                  2*p1k1*sqrtS**2*(4*p1k2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) - 
     &                     sqrtS**2*(-4*muu4**2*x1 + muu2**2*(4 + 4*muu3**2 + 4*muu4**2 - 4*x1 - 6*x2) - 4*muu1**2*(muu4**2 - x2) + x2 + muu3**2*x2 + muu4**2*x2 + 3*x1*x2)) + 
     &                  sqrtS**4*(-((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2) + 2*p1k2*(1 + 4*muu1**2 - 2*muu2**2 + muu4**2 - x1)*x1 - 
     &                     sqrtS**2*(-((-1 + muu2**2)*(-1 + 2*muu2**2 - muu4**2 + x1)*(1 + muu4**2 - x1 - x2)) + 
     &                        muu1**2*(1 + muu4**4 - 2*x1 + x1**2 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2 + 5*x1*x2 - muu4**2*(-2 + 6*x1 + x2))) + 
     &                     muu3**2*(p1k2*(8 - 8*muu1**2 + 2*x1) + sqrtS**2*(muu1**2*(-2 - 2*muu2**2 - 6*muu4**2 + 2*x1 + 5*x2) + (-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + 2*x1 + 5*x2)))))
     &                 *HL(indj,i,tt,gen) - 2*muu2*muu3*sqrtS**2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HR(indj,i,tt,gen)))) + 
     &      dconjg(ER(indi,j,tt,gen))*(-2*sqrtS**2*dconjg(FL(indj,j,tt,gen))*
     &          (muu1*GL(indi,i,tt,gen)*(muu3*(-4*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + 
     &                  sqrtS**2*(2*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) - 
     &                     sqrtS**2*(2*muu2**2*(1 + muu3**2 - x1 - x2) + (-1 + 2*muu1**2 + muu3**2 + x1)*x2 + muu4**2*(2 - 2*muu1**2 - 2*x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 - 2*p1k2*x1 - sqrtS**2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 
     &                  sqrtS**2*x1**2 + muu3**2*(4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - 2*x1 - 3*x2)) + muu4**2*(-4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 2*x1 - x2)) + 
     &                  2*p1k1*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2)*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*(muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu2**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 - 2*p1k2*x1 + 2*sqrtS**2*x1 + 
     &                  muu1**2*(4*p1k2 + sqrtS**2*(2 + 2*muu3**2 + 2*muu4**2 - 2*x1 - 3*x2)) + muu2**2*(-4*p1k1 + sqrtS**2*(4 + 2*muu1**2 + 2*muu3**2 + 2*muu4**2 - 2*x1 - x2)) + 
     &                  2*p1k1*x2 + sqrtS**2*x2)*HL(indj,i,tt,gen) + 2*muu2*(4*p1k1**2 + muu3**2*(-2*p1k2*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + muu1**2*sqrtS**4*(1 + muu4**2 - x2) + 
     &                  p1k1*sqrtS**2*(-3 - 2*muu1**2 + 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + x1 + 2*x2))*HR(indj,i,tt,gen))) + 
     &         dconjg(FR(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             (-2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*HL(indj,i,tt,gen) + 
     &               (-2*p1k1*(4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + 
     &                  sqrtS**2*(2*p1k2*x1 + muu3**2*(4*p1k2 + sqrtS**2*(1 - muu1**2 + 3*muu2**2 + x1 - 2*x2)) - 
     &                     sqrtS**2*(-1 - muu4**2 - muu4**2*x1 + x1**2 + muu1**2*(1 + muu4**2 - x1 - x2) + muu2**2*(1 + muu4**2 + x1 - x2) + x2)))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(-2*muu2*muu3*sqrtS**2*(p1k1*(8*p1k2 + 2*sqrtS**2*(2*muu4**2 - x2)) + 
     &                  sqrtS**2*(-2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*p1k2 + 
     &                     sqrtS**2*((-1 + muu2**2)*(1 + muu3**2 + muu4**2 - x1 - x2) + muu1**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &               (16*p1k1**2*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) - 2*p1k1*sqrtS**2*
     &                   (4*p1k2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x1) - 
     &                     sqrtS**2*(-4*muu4**2*x1 + muu2**2*(4 + 4*muu3**2 + 4*muu4**2 - 4*x1 - 6*x2) - 4*muu1**2*(muu4**2 - x2) + x2 + muu3**2*x2 + muu4**2*x2 + 3*x1*x2)) + 
     &                  sqrtS**4*(-((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2) + 2*p1k2*(1 + 4*muu1**2 - 2*muu2**2 + muu4**2 - x1)*x1 - 
     &                     sqrtS**2*(-((-1 + muu2**2)*(-1 + 2*muu2**2 - muu4**2 + x1)*(1 + muu4**2 - x1 - x2)) + 
     &                        muu1**2*(1 + muu4**4 - 2*x1 + x1**2 + 2*muu2**2*(1 + muu4**2 - x1 - x2) - x2 + 5*x1*x2 - muu4**2*(-2 + 6*x1 + x2))) + 
     &                     muu3**2*(p1k2*(8 - 8*muu1**2 + 2*x1) + sqrtS**2*(muu1**2*(-2 - 2*muu2**2 - 6*muu4**2 + 2*x1 + 5*x2) + (-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + 2*x1 + 5*x2)))))
     &                 *HR(indj,i,tt,gen))))))/
     &  (2.*sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))  )

** TUchannel bremsstrahlung 31

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(   (CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*(-2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*
     &             ((2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + 2*muu2**2*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     x1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + muu4**2 + x2))))*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 + 2*muu3**2 + 2*muu4**2 - x1 - x2))*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((-2*p1k1*(16*p1k2**2 - 4*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2) + 
     &                     sqrtS**4*((1 + 4*muu1**2 - 2*muu2**2 + muu3**2 - x2)*x2 + muu4**2*(4 - 4*muu1**2 + x2))) + 
     &                  sqrtS**2*(-sqrtS**4 + muu1**2*sqrtS**4 + (-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + 2*muu1**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + 
     &                     muu1**2*muu4**4*sqrtS**4 + 16*p1k2**2*x1 - 2*p1k2*sqrtS**2*x1 - 8*muu1**2*p1k2*sqrtS**2*x1 - 2*muu4**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - muu1**2*sqrtS**4*x1 + 
     &                     5*muu4**2*sqrtS**4*x1 - 5*muu1**2*muu4**2*sqrtS**4*x1 - 2*muu2**4*sqrtS**4*(1 + muu4**2 - x1 - x2) + 2*sqrtS**4*x2 - 2*muu1**2*sqrtS**4*x2 + 
     &                     2*muu4**2*sqrtS**4*x2 - 2*muu1**2*muu4**2*sqrtS**4*x2 - 6*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 + 5*muu1**2*sqrtS**4*x1*x2 - sqrtS**4*x2**2 + 
     &                     muu1**2*sqrtS**4*x2**2 + muu2**2*sqrtS**2*(p1k2*(-8 - 8*muu4**2 + 12*x1 + 8*x2) + 
     &                        sqrtS**2*(muu4**4 + muu4**2*(4 - 5*x1 - 2*x2) + 2*muu1**2*(1 + muu4**2 - x1 - x2) + (-3 + x2)*(-1 + x1 + x2))) - 
     &                     muu3**2*(16*p1k2**2 + 2*p1k2*sqrtS**2*(-4*muu1**2 + 4*muu2**2 + x1 - 4*x2) + 
     &                        sqrtS**4*((-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + x1 + 2*x2) + muu1**2*(-2 - 2*muu2**2 - 6*muu4**2 + x1 + 6*x2)))))*HL(indj,i,tt,gen) + 
     &               4*muu2*muu3*sqrtS**2*(4*p1k2**2 + muu4**2*(-2*p1k1*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + muu1**2*sqrtS**4*(1 + muu3**2 - x1) + 
     &                  p1k2*sqrtS**2*(-3 - 2*muu1**2 + 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + 2*x1 + x2))*HR(indj,i,tt,gen))) + 
     &         2*sqrtS**2*dconjg(FR(indj,j,tt,gen))*(muu4*GL(indi,i,tt,gen)*(muu2*
     &                (-2*p1k1*(-4*p1k2 + (1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*sqrtS**2) + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                HL(indj,i,tt,gen) + muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu2**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 + 2*p1k2*x1 + sqrtS**2*x1 + 
     &                  muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) + muu2**2*(-4*p1k2 + sqrtS**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) - 
     &                  2*p1k1*x2 + 2*sqrtS**2*x2)*HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*
     &             (muu2*sqrtS**2*(muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 
     &                  muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - 
     &                  sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + sqrtS**2*x2**2)*HL(indj,i,tt,gen) + 
     &               muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + sqrtS**2*
     &                   (2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + sqrtS**2*
     &                      (-1 - muu3**2 - muu4**2 + x1 + 2*muu4**2*x1 + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - muu3**2*x2 - muu4**2*x2 + x2**2 + 
     &                        muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HR(indj,i,tt,gen)))) + 
     &      dconjg(ER(indi,j,tt,gen))*(2*sqrtS**2*dconjg(FL(indj,j,tt,gen))*(muu4*GR(indi,i,tt,gen)*
     &             (muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu2**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 + 2*p1k2*x1 + sqrtS**2*x1 + 
     &                  muu1**2*(4*p1k1 + sqrtS**2*(2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 3*x1 - 2*x2)) + muu2**2*(-4*p1k2 + sqrtS**2*(4 + 2*muu3**2 + 2*muu4**2 - x1 - 2*x2)) - 
     &                  2*p1k1*x2 + 2*sqrtS**2*x2)*HL(indj,i,tt,gen) + muu2*(-2*p1k1*(-4*p1k2 + (1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2)*sqrtS**2) + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))))*
     &                HR(indj,i,tt,gen)) + muu1*GL(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                  sqrtS**2*(2*p1k2*(1 + muu3**2 + muu4**2 - x1 - x2) + 
     &                     sqrtS**2*(-1 - muu3**2 - muu4**2 + x1 + 2*muu4**2*x1 + muu1**2*(1 + muu3**2 + muu4**2 - x1 - x2) - muu3**2*x2 - muu4**2*x2 + x2**2 + 
     &                        muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &               muu2*sqrtS**2*(muu1**2*sqrtS**2 - muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - muu1**2*sqrtS**2*x1 + muu2**2*sqrtS**2*x1 + 
     &                  muu4**2*(4*p1k1 + sqrtS**2*(1 + muu1**2 - muu2**2 + 6*muu3**2 - 3*x1 - 2*x2)) + muu3**2*(-4*p1k2 + sqrtS**2*(1 + muu1**2 - muu2**2 - x1 - 2*x2)) - 2*p1k1*x2 - 
     &                  sqrtS**2*x2 - muu1**2*sqrtS**2*x2 + muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + sqrtS**2*x2**2)*HR(indj,i,tt,gen))) + 
     &         dconjg(FR(indj,j,tt,gen))*(-2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &             (-2*muu2*muu3*sqrtS**2*(-2*p1k1 + 2*p1k2 + sqrtS**2*(-1 - muu1**2 + muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2))*HL(indj,i,tt,gen) + 
     &               (2*p1k1*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2)) + 
     &                  sqrtS**2*(muu3**2*(4*p1k2 + sqrtS**2*(2 - 2*muu1**2 + x1 - 2*x2)) + 2*muu2**2*sqrtS**2*(1 + muu4**2 - x1 - x2) + 
     &                     x1*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + muu4**2 + x2))))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(4*muu2*muu3*sqrtS**2*(4*p1k2**2 + muu4**2*(-2*p1k1*sqrtS**2 + sqrtS**4 - muu2**2*sqrtS**4) + muu1**2*sqrtS**4*(1 + muu3**2 - x1) + 
     &                  p1k2*sqrtS**2*(-3 - 2*muu1**2 + 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + 2*x1 + x2))*HL(indj,i,tt,gen) + 
     &               (-2*p1k1*(16*p1k2**2 - 4*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 + x2) + 
     &                     sqrtS**4*((1 + 4*muu1**2 - 2*muu2**2 + muu3**2 - x2)*x2 + muu4**2*(4 - 4*muu1**2 + x2))) + 
     &                  sqrtS**2*(-sqrtS**4 + muu1**2*sqrtS**4 + (-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + 2*muu1**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + 
     &                     muu1**2*muu4**4*sqrtS**4 + 16*p1k2**2*x1 - 2*p1k2*sqrtS**2*x1 - 8*muu1**2*p1k2*sqrtS**2*x1 - 2*muu4**2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - muu1**2*sqrtS**4*x1 + 
     &                     5*muu4**2*sqrtS**4*x1 - 5*muu1**2*muu4**2*sqrtS**4*x1 - 2*muu2**4*sqrtS**4*(1 + muu4**2 - x1 - x2) + 2*sqrtS**4*x2 - 2*muu1**2*sqrtS**4*x2 + 
     &                     2*muu4**2*sqrtS**4*x2 - 2*muu1**2*muu4**2*sqrtS**4*x2 - 6*p1k2*sqrtS**2*x1*x2 - sqrtS**4*x1*x2 + 5*muu1**2*sqrtS**4*x1*x2 - sqrtS**4*x2**2 + 
     &                     muu1**2*sqrtS**4*x2**2 + muu2**2*sqrtS**2*(p1k2*(-8 - 8*muu4**2 + 12*x1 + 8*x2) + 
     &                        sqrtS**2*(muu4**4 + muu4**2*(4 - 5*x1 - 2*x2) + 2*muu1**2*(1 + muu4**2 - x1 - x2) + (-3 + x2)*(-1 + x1 + x2))) - 
     &                     muu3**2*(16*p1k2**2 + 2*p1k2*sqrtS**2*(-4*muu1**2 + 4*muu2**2 + x1 - 4*x2) + 
     &                        sqrtS**4*((-1 + muu2**2)*(-2 + 2*muu2**2 - 6*muu4**2 + x1 + 2*x2) + muu1**2*(-2 - 2*muu2**2 - 6*muu4**2 + x1 + 6*x2)))))*HR(indj,i,tt,gen))))))/
     &  (2.*sqrtS**2*(-1 + muu3**2 - muu4**2 + x2)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))   )

** TUchannel bremsstrahlung 32

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(  -(CF*gs**2*NC*(dconjg(ELc(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     &           (-2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*((-8*p1k1**2 + sqrtS**4*(x1*(3 - 2*muu1**2 + 3*muu4**2 - 3*x1 - 2*x2) - 2*muu2**2*(1 + muu4**2 - x1 - x2)) + 
     &                   2*p1k1*sqrtS**2*(-3 + 2*muu1**2 - 2*muu2**2 - 3*muu3**2 - 3*muu4**2 + 5*x1 + 2*x2) + muu3**2*(-4*p1k2*sqrtS**2 + sqrtS**4*(-2 + 2*muu1**2 + 3*x1 + 2*x2)))*
     &                 HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*HR(indj,i,tt,gen)) + 
     &             GL(indi,i,tt,gen)*((16*p1k1**2*(2*p1k2 - muu4**2*sqrtS**2) + 
     &                   sqrtS**4*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-1 + 2*muu1**2 - muu4**2 + x1) + 
     &                      muu3**2*(-2*p1k2*(-4 + 4*muu2**2 + x1) + sqrtS**2*
     &                          (-2*muu1**4 + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - 6*x1 - 5*x2) + (-1 + muu2**2)*(2 + 6*muu4**2 - 2*x1 - x2))) - 
     &                      sqrtS**2*(2*muu1**4 - muu1**2*(3 + 2*muu2**2 + muu4**2 - 5*x1) - (-1 + muu2**2)*(1 + muu4**2 - x1))*(1 + muu4**2 - x1 - x2)) - 
     &                   2*p1k1*sqrtS**2*(4*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 - muu3**2 - muu4**2 + 3*x1) + 
     &                      sqrtS**2*(4*muu2**2*muu4**2 - 4*muu4**2*x1 + x2 + muu3**2*x2 + muu4**2*x2 - x1*x2 + muu1**2*(-4 - 4*muu3**2 - 4*muu4**2 + 4*x1 + 2*x2))))*HL(indj,i,tt,gen) + 
     &                4*muu2*muu3*sqrtS**2*(p1k1*(-4*p1k2 + 2*muu4**2*sqrtS**2) + 
     &                   sqrtS**2*(-(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1))) + p1k2*(1 + 2*muu1**2 - 2*muu2**2 - muu3**2 - muu4**2 + x1)))*HR(indj,i,tt,gen))) + 
     &          2*sqrtS**2*dconjg(FR(indj,j,tt,gen))*(muu4*GL(indi,i,tt,gen)*
     &              (muu2*(sqrtS**2*(muu3**2*(4*p1k2 + (-1 - 3*muu1**2 + muu2**2)*sqrtS**2) - 2*p1k2*x1 + 
     &                      sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) + 2*p1k1*(4*p1k2 - sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*
     &                 HL(indj,i,tt,gen) + muu3*sqrtS**2*(-2*sqrtS**2 - 2*muu1**4*sqrtS**2 - 2*muu3**2*sqrtS**2 - 2*muu4**2*sqrtS**2 + 2*p1k2*x1 + 2*sqrtS**2*x1 + 
     &                   muu1**2*(4*p1k1 + sqrtS**2*(4 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 4*x1 - x2)) + muu2**2*(-4*p1k2 + sqrtS**2*(2 + 2*muu3**2 + 2*muu4**2 - 2*x1 - x2)) - 2*p1k1*x2 + 
     &                   sqrtS**2*x2)*HR(indj,i,tt,gen)) + muu1*GR(indi,i,tt,gen)*
     &              (muu2*sqrtS**2*(-(muu1**2*sqrtS**2) + muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 
     &                   sqrtS**2*x1**2 - 2*p1k1*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + 
     &                   muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HL(indj,i,tt,gen) + 
     &                muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + sqrtS**2*
     &                    (2*p1k2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1) + 
     &                      sqrtS**2*(-1 - muu3**2 - muu4**2 + x1 + 2*muu4**2*x1 + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - muu3**2*x2 - muu4**2*x2 + x1*x2 + 
     &                         muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HR(indj,i,tt,gen)))) + 
     &       dconjg(ER(indi,j,tt,gen))*(-2*sqrtS**2*dconjg(FL(indj,j,tt,gen))*
     &           (muu4*GR(indi,i,tt,gen)*(muu3*sqrtS**2*(2*sqrtS**2 + 2*muu1**4*sqrtS**2 + 2*muu3**2*sqrtS**2 + 2*muu4**2*sqrtS**2 - 2*p1k2*x1 - 2*sqrtS**2*x1 - 
     &                   muu1**2*(4*p1k1 + sqrtS**2*(4 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - 4*x1 - x2)) + 2*p1k1*x2 - sqrtS**2*x2 + 
     &                   muu2**2*(4*p1k2 + sqrtS**2*(-2 - 2*muu3**2 - 2*muu4**2 + 2*x1 + x2)))*HL(indj,i,tt,gen) + 
     &                muu2*(sqrtS**2*(muu3**2*(-4*p1k2 + (1 + 3*muu1**2 - muu2**2)*sqrtS**2) + 2*p1k2*x1 - 
     &                      sqrtS**2*((-1 + muu2**2)*(1 + muu4**2 - x1 - x2) + muu1**2*(1 + muu4**2 + x1 - x2))) + 2*p1k1*(-4*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - x1 + x2)))*
     &                 HR(indj,i,tt,gen)) - muu1*GL(indi,i,tt,gen)*(muu3*(p1k1*(8*p1k2 - 2*sqrtS**2*(2*muu4**2 + x2)) + 
     &                   sqrtS**2*(2*p1k2*(1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - 2*x1) + 
     &                      sqrtS**2*(-1 - muu3**2 - muu4**2 + x1 + 2*muu4**2*x1 + muu2**2*(1 + muu3**2 - 3*muu4**2 - x1 - x2) - muu3**2*x2 - muu4**2*x2 + x1*x2 + 
     &                         muu1**2*(1 + muu3**2 + muu4**2 - x1 + x2))))*HL(indj,i,tt,gen) + 
     &                muu2*sqrtS**2*(-(muu1**2*sqrtS**2) + muu2**2*sqrtS**2 + muu3**4*sqrtS**2 + muu4**4*sqrtS**2 + 2*p1k2*x1 - sqrtS**2*x1 + muu1**2*sqrtS**2*x1 - muu2**2*sqrtS**2*x1 + 
     &                   sqrtS**2*x1**2 - 2*p1k1*x2 + muu1**2*sqrtS**2*x2 - muu2**2*sqrtS**2*x2 + sqrtS**2*x1*x2 + muu3**2*(-4*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + 2*x1 + x2)) + 
     &                   muu4**2*(4*p1k1 - sqrtS**2*(-1 + muu1**2 - muu2**2 - 6*muu3**2 + 4*x1 + x2)))*HR(indj,i,tt,gen))) + 
     &          dconjg(FR(indj,j,tt,gen))*(2*muu1*muu4*sqrtS**2*GL(indi,i,tt,gen)*
     &              (-2*muu2*muu3*sqrtS**2*(2*p1k1 - 2*p1k2 + sqrtS**2*(1 - muu1**2 + muu2**2 + 2*muu3**2 + 2*muu4**2 - 2*x1))*HL(indj,i,tt,gen) + 
     &                (8*p1k1**2 + 2*p1k1*sqrtS**2*(3 - 2*muu1**2 + 2*muu2**2 + 3*muu3**2 + 3*muu4**2 - 5*x1 - 2*x2) + muu3**2*(4*p1k2*sqrtS**2 - sqrtS**4*(-2 + 2*muu1**2 + 3*x1 + 2*x2)) + 
     &                   sqrtS**4*(2*muu2**2*(1 + muu4**2 - x1 - x2) + x1*(-3 + 2*muu1**2 - 3*muu4**2 + 3*x1 + 2*x2)))*HR(indj,i,tt,gen)) + 
     &             GR(indi,i,tt,gen)*(-4*muu2*muu3*sqrtS**2*(p1k1*(4*p1k2 - 2*muu4**2*sqrtS**2) + 
     &                   sqrtS**2*(sqrtS**2*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 + muu3**2 - x1)) + p1k2*(-1 - 2*muu1**2 + 2*muu2**2 + muu3**2 + muu4**2 - x1)))*HL(indj,i,tt,gen) + 
     &                (16*p1k1**2*(2*p1k2 - muu4**2*sqrtS**2) + sqrtS**4*((-1 + muu1**2 + muu2**2)*muu3**4*sqrtS**2 + 2*p1k2*x1*(-1 + 2*muu1**2 - muu4**2 + x1) + 
     &                      muu3**2*(-2*p1k2*(-4 + 4*muu2**2 + x1) + sqrtS**2*
     &                          (-2*muu1**4 + muu1**2*(4 + 2*muu2**2 + 6*muu4**2 - 6*x1 - 5*x2) + (-1 + muu2**2)*(2 + 6*muu4**2 - 2*x1 - x2))) - 
     &                      sqrtS**2*(2*muu1**4 - muu1**2*(3 + 2*muu2**2 + muu4**2 - 5*x1) - (-1 + muu2**2)*(1 + muu4**2 - x1))*(1 + muu4**2 - x1 - x2)) - 
     &                   2*p1k1*sqrtS**2*(4*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 - muu3**2 - muu4**2 + 3*x1) + 
     &                      sqrtS**2*(4*muu2**2*muu4**2 - 4*muu4**2*x1 + x2 + muu3**2*x2 + muu4**2*x2 - x1*x2 + muu1**2*(-4 - 4*muu3**2 - 4*muu4**2 + 4*x1 + 2*x2))))*HR(indj,i,tt,gen))))))/
     &  (2.*sqrtS**2*(1 + muu3**2 - muu4**2 - x1)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (-2*p1k1 + (muu1**2 + muu3**2)*sqrtS**2 - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))    )

** TUchannel bremsstrahlung 33

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU + 2d0*dreal(  (CF*gs**2*NC*sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + x1 + x2)*
     &    (dconjg(ER(indi,j,tt,gen))*(dconjg(FR(indj,j,tt,gen))*(-4*muu2*muu3*sqrtS**2*(muu1*muu4*sqrtS**2*GL(indi,i,tt,gen) + p1k2*GR(indi,i,tt,gen))*HL(indj,i,tt,gen) + 
     &            (-2*muu1*muu4*sqrtS**2*(-2*p1k1 + sqrtS**2*x1)*GL(indi,i,tt,gen) + 
     &               (sqrtS**2*(-2*p1k2*x1 + (-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + p1k1*(8*p1k2 - 2*sqrtS**2*x2))*GR(indi,i,tt,gen))*HR(indj,i,tt,gen)) + 
     &         2*sqrtS**2*dconjg(FL(indj,j,tt,gen))*(muu4*GR(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HL(indj,i,tt,gen) - 2*muu2*p1k1*HR(indj,i,tt,gen)) + 
     &            muu1*GL(indi,i,tt,gen)*(muu3*(2*p1k2 - sqrtS**2*x2)*HL(indj,i,tt,gen) + muu2*sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)*HR(indj,i,tt,gen)))) + 
     &      dconjg(ELc(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*(-2*muu1*muu4*sqrtS**2*GR(indi,i,tt,gen)*
     &             ((-2*p1k1 + sqrtS**2*x1)*HL(indj,i,tt,gen) + 2*muu2*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((sqrtS**2*(-2*p1k2*x1 + (-1 + muu1**2 + muu2**2)*sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)) + p1k1*(8*p1k2 - 2*sqrtS**2*x2))*HL(indj,i,tt,gen) - 
     &               4*muu2*muu3*p1k2*sqrtS**2*HR(indj,i,tt,gen))) + 2*sqrtS**2*dconjg(FR(indj,j,tt,gen))*
     &          (muu4*GL(indi,i,tt,gen)*(-2*muu2*p1k1*HL(indj,i,tt,gen) + (-1 + muu1**2 + muu2**2)*muu3*sqrtS**2*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(muu2*sqrtS**2*(1 + muu3**2 + muu4**2 - x1 - x2)*HL(indj,i,tt,gen) + muu3*(2*p1k2 - sqrtS**2*x2)*HR(indj,i,tt,gen))))))/
     &  (2.*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2 - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    (2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen)))*(-2*p1k2 - sqrtS**2*(muu2**2 + muu4**2 - x2) + (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))   )

        if (cutoff.lt.zerocheck) then

** TUchannel dipole 312

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU - 2d0*dreal(  (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     &         (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))*
     &    (dconjg(ELc(indi,j,tt,gen))*(dconjg(FR(indj,j,tt,gen))*(2*muu1*muu2*(-1 + muu3**2 + muu4**2)*sqrtS**6*(4*muu4**2 - x2**2)*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &            muu1*(-1 + muu1**2 - muu2**2)*muu3*(1 - muu3**2 + muu4**2)*sqrtS**6*(4*muu4**2 - x2**2)*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            4*muu1*muu3*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            muu1*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            muu1**3*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu1*muu2**2*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*
     &             HR(indj,i,tt,gen) + muu4*GL(indi,i,tt,gen)*(muu2*(-4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 + (1 + muu3**2)*sqrtS**4*x2**2 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &                     Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))))*HL(indj,i,tt,gen) + 
     &               2*(-1 + muu1**2 + muu2**2)*muu3*sqrtS**6*(4*muu4**2 - x2**2)*HR(indj,i,tt,gen))) + 
     &         dconjg(FL(indj,j,tt,gen))*(muu1*muu4*GR(indi,i,tt,gen)*((4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 + 
     &                  (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) + 
     &                  4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                  (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HL(indj,i,tt,gen)
     &                + 4*muu2*muu3*sqrtS**6*(-4*muu4**2 + x2**2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((8*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*p1k2**2*sqrtS**2 - 
     &                  4*p1k2*((1 + muu1**2 - muu2**2)*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4*x2 + 
     &                     (muu1 - muu2)*(muu1 + muu2)*(muu3 - muu4)*(muu3 + muu4)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + sqrtS**2*
     &                   (-2*(-1 + (muu1**2 - muu2**2)**2)*muu4**6*sqrtS**4 + 
     &                     muu4**4*sqrtS**4*(4*(-1 + muu1**2 + muu2**2 + (-1 + (muu1**2 - muu2**2)**2)*muu3**2) + (muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*x2**2) + 
     &                     x2*(muu1**4*muu3**2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + 
     &                        muu2**2*(-1 + muu2**2)*muu3**2*((-1 + muu3**2)*sqrtS**4*x2 + 
     &                           Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + 
     &                        muu1**2*(-((-1 + muu3**2)*(2 + (-1 + 2*muu2**2)*muu3**2)*sqrtS**4*x2) + 
     &                           (1 - 2*muu2**2)*muu3**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))) - 
     &                     muu4**2*((muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*
     &                         Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + sqrtS**4*
     &                         (muu1**4*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2)) + 
     &                           muu1**2*(4 + 3*x2**2 + 2*muu3**2*(-2 + x2**2) - 2*muu2**2*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2))) + 
     &                           (-1 + muu2**2)*(2*(-1 + muu3**2)**2 + muu2**2*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2)))))))*HL(indj,i,tt,gen) + 
     &               muu2*muu3*(4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                     x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HR(indj,i,tt,gen))
     &            )) + dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     &          (2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**6*(4*muu4**2 - x2**2)*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &            muu2*(1 + muu1**2 - muu2**2)*muu4*(-1 - muu3**2 + muu4**2)*sqrtS**6*(4*muu4**2 - x2**2)*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            4*muu2*muu4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            muu1**2*muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*GR(indi,i,tt,gen)*
     &             HR(indj,i,tt,gen) - muu2**3*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))*
     &             GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + muu1*GL(indi,i,tt,gen)*
     &             (muu3*(4*(-1 + muu1**2 - muu2**2)*muu4**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu3**2)*sqrtS**6*x2**2 - 
     &                  (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(-4 + 4*muu3**2 + x2**2) - 
     &                  4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HL(indj,i,tt,gen) + 
     &               2*muu2*(-1 + muu3**2 + muu4**2)*sqrtS**6*(4*muu4**2 - x2**2)*HR(indj,i,tt,gen))) + 
     &         dconjg(FR(indj,j,tt,gen))*(muu1*muu4*GL(indi,i,tt,gen)*(4*muu2*muu3*sqrtS**6*(-4*muu4**2 + x2**2)*HL(indj,i,tt,gen) + 
     &               (4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 + 
     &                  (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) + 
     &                  4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &                  (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))*HR(indj,i,tt,gen))
     &             + GR(indi,i,tt,gen)*(muu2*muu3*(4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(-4 + 4*muu3**2 + x2**2) + 
     &                     x2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))*HL(indj,i,tt,gen)
     &                + (8*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*p1k2**2*sqrtS**2 - 
     &                  4*p1k2*((1 + muu1**2 - muu2**2)*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4*x2 + 
     &                     (muu1 - muu2)*(muu1 + muu2)*(muu3 - muu4)*(muu3 + muu4)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + sqrtS**2*
     &                   (-2*(-1 + (muu1**2 - muu2**2)**2)*muu4**6*sqrtS**4 + 
     &                     muu4**4*sqrtS**4*(4*(-1 + muu1**2 + muu2**2 + (-1 + (muu1**2 - muu2**2)**2)*muu3**2) + (muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*x2**2) + 
     &                     x2*(muu1**4*muu3**2*((-1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + 
     &                        muu2**2*(-1 + muu2**2)*muu3**2*((-1 + muu3**2)*sqrtS**4*x2 + 
     &                           Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2))) + 
     &                        muu1**2*(-((-1 + muu3**2)*(2 + (-1 + 2*muu2**2)*muu3**2)*sqrtS**4*x2) + 
     &                           (1 - 2*muu2**2)*muu3**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))) - 
     &                     muu4**2*((muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*
     &                         Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + sqrtS**4*
     &                         (muu1**4*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2)) + 
     &                           muu1**2*(4 + 3*x2**2 + 2*muu3**2*(-2 + x2**2) - 2*muu2**2*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2))) + 
     &                           (-1 + muu2**2)*(2*(-1 + muu3**2)**2 + muu2**2*(-2 + x2**2 + 2*muu3**2*(muu3**2 + x2**2)))))))*HR(indj,i,tt,gen))))))/
     &  (sqrtS**4*(1 - muu3**2 + muu4**2 - x2)*(4*muu4**2 - x2**2)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (4*(1 - muu1**2 + muu2**2)*muu4**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(1 + muu3**2)*sqrtS**6*x2**2 + (-1 + muu1**2 - muu2**2)*muu4**2*sqrtS**6*(4 + 4*muu3**2 + x2**2) + 
     &         4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) + 
     &         (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x2*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))/
     &       (2.*sqrtS**4*(4*muu4**2 - x2**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*((muu1**2 + muu3**2)*sqrtS**2 - 
     &      (4*p1k2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)) - 
     &         (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu4**4*sqrtS**4 - muu4**2*sqrtS**4*(4 + 4*muu3**2 + x2**2) + 
     &            x2*((1 + muu3**2)*sqrtS**4*x2 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu4**2 + x2**2)))))/
     &       (2.*sqrtS**4*(4*muu4**2 - x2**2)) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))))  )

** TUchannel dipole 321 

	ChiChi2QQ_RadiationTU = ChiChi2QQ_RadiationTU - 2d0*dreal(  (CF*GS**2*NC*((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) + (Sqrt((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*
     &         (-3*muu3**4 + muu4**4 + muu4**2*(x1 - x2) - (-1 + x1)*(-1 + 2*x1 + x2) + muu3**2*(-4 - 2*muu4**2 + 5*x1 + x2)))/
     &       ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))*
     &    (dconjg(ELc(indi,j,tt,gen))*(dconjg(FR(indj,j,tt,gen))*(-(muu2*(1 + muu1**2 - muu2**2)*muu4*(1 + muu3**2 - muu4**2)*sqrtS**6*(4*muu3**2 - x1**2)*GL(indi,i,tt,gen)*
     &               HL(indj,i,tt,gen)) + 4*muu2*muu4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GL(indi,i,tt,gen)*
     &             HL(indj,i,tt,gen) - muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GL(indi,i,tt,gen)*
     &             HL(indj,i,tt,gen) - muu1**2*muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*
     &             GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + muu2**3*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*
     &             Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GL(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &            2*muu1*muu2*(-1 + muu3**2 + muu4**2)*sqrtS**6*(4*muu3**2 - x1**2)*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) + 
     &            2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**6*(4*muu3**2 - x1**2)*GL(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu1*(-1 + muu1**2 - muu2**2)*muu3*(-1 + muu3**2 - muu4**2)*sqrtS**6*(4*muu3**2 - x1**2)*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            4*muu1*muu3*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu1*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu1**3*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            muu1*muu2**2*muu3*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*
     &             HR(indj,i,tt,gen)) + dconjg(FL(indj,j,tt,gen))*(muu1*muu4*GR(indi,i,tt,gen)*
     &             ((4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &                  4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                  (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*HL(indj,i,tt,gen) + 4*muu2*muu3*sqrtS**6*(-4*muu3**2 + x1**2)*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*((8*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*p1k1**2*sqrtS**2 - 
     &                  4*p1k1*((1 + muu1**2 - muu2**2)*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4*x1 - 
     &                     (muu1 - muu2)*(muu1 + muu2)*(muu3 - muu4)*(muu3 + muu4)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + sqrtS**2*
     &                   (-2*(-1 + (muu1**2 - muu2**2)**2)*muu3**6*sqrtS**4 + 
     &                     muu3**4*sqrtS**4*(4*(-1 + muu1**2 + muu2**2 + (-1 + (muu1**2 - muu2**2)**2)*muu4**2) + (muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*x1**2) - 
     &                     muu3**2*((muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*
     &                         Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + sqrtS**4*
     &                         (2*(-1 + muu4**2)*((-1 + muu2**2)**2 + (-1 + muu2**4)*muu4**2 + muu1**4*(1 + muu4**2) - 2*muu1**2*(1 + muu2**2*(1 + muu4**2))) + 
     &                           (3*muu1**2 + muu1**4 - muu2**2 - 2*muu1**2*muu2**2 + muu2**4 + 2*(muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*muu4**2)*x1**2)) + 
     &                     x1*(muu1**4*muu4**2*((-1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &                        muu2**2*(-1 + muu2**2)*muu4**2*((-1 + muu4**2)*sqrtS**4*x1 + 
     &                           Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &                        muu1**2*(-((-1 + muu4**2)*(2 + (-1 + 2*muu2**2)*muu4**2)*sqrtS**4*x1) + 
     &                           (1 - 2*muu2**2)*muu4**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))))*HL(indj,i,tt,gen) + 
     &               muu2*muu3*(-4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                     x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*HR(indj,i,tt,gen)))
     &         ) + dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*(2*(-1 + muu1**2 + muu2**2)*muu3*muu4*sqrtS**6*(4*muu3**2 - x1**2)*GR(indi,i,tt,gen)*HL(indj,i,tt,gen) - 
     &            muu2*(1 + muu1**2 - muu2**2)*muu4*(1 + muu3**2 - muu4**2)*sqrtS**6*(4*muu3**2 - x1**2)*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + 
     &            4*muu2*muu4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*HR(indj,i,tt,gen) - 
     &            muu1**2*muu2*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*GR(indi,i,tt,gen)*
     &             HR(indj,i,tt,gen) + muu2**3*muu4*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))*
     &             GR(indi,i,tt,gen)*HR(indj,i,tt,gen) + muu1*GL(indi,i,tt,gen)*
     &             (muu3*(4*(1 - muu1**2 + muu2**2)*muu3**4*sqrtS**6 - (-1 + muu1**2 - muu2**2)*(1 + muu4**2)*sqrtS**6*x1**2 + 
     &                  4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (-1 - muu1**2 + muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(4 + 4*muu4**2 + x1**2))*HL(indj,i,tt,gen) + 
     &               2*muu2*(-1 + muu3**2 + muu4**2)*sqrtS**6*(4*muu3**2 - x1**2)*HR(indj,i,tt,gen))) + 
     &         dconjg(FR(indj,j,tt,gen))*(muu1*muu4*GL(indi,i,tt,gen)*(4*muu2*muu3*sqrtS**6*(-4*muu3**2 + x1**2)*HL(indj,i,tt,gen) + 
     &               (4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &                  4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &                  (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(muu2*muu3*(-4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &                  (1 + muu1**2 - muu2**2)*sqrtS**2*(4*muu3**4*sqrtS**4 - muu3**2*sqrtS**4*(4 + 4*muu4**2 + x1**2) + 
     &                     x1*((1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))))*HL(indj,i,tt,gen)
     &                + (8*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*p1k1**2*sqrtS**2 - 
     &                  4*p1k1*((1 + muu1**2 - muu2**2)*((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4*x1 - 
     &                     (muu1 - muu2)*(muu1 + muu2)*(muu3 - muu4)*(muu3 + muu4)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*
     &                      Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + sqrtS**2*
     &                   (-2*(-1 + (muu1**2 - muu2**2)**2)*muu3**6*sqrtS**4 + 
     &                     muu3**4*sqrtS**4*(4*(-1 + muu1**2 + muu2**2 + (-1 + (muu1**2 - muu2**2)**2)*muu4**2) + (muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*x1**2) - 
     &                     muu3**2*((muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*
     &                         Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + sqrtS**4*
     &                         (2*(-1 + muu4**2)*((-1 + muu2**2)**2 + (-1 + muu2**4)*muu4**2 + muu1**4*(1 + muu4**2) - 2*muu1**2*(1 + muu2**2*(1 + muu4**2))) + 
     &                           (3*muu1**2 + muu1**4 - muu2**2 - 2*muu1**2*muu2**2 + muu2**4 + 2*(muu1 - muu2)*(muu1 + muu2)*(1 + muu1**2 - muu2**2)*muu4**2)*x1**2)) + 
     &                     x1*(muu1**4*muu4**2*((-1 + muu4**2)*sqrtS**4*x1 + Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &                        muu2**2*(-1 + muu2**2)*muu4**2*((-1 + muu4**2)*sqrtS**4*x1 + 
     &                           Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))) + 
     &                        muu1**2*(-((-1 + muu4**2)*(2 + (-1 + 2*muu2**2)*muu4**2)*sqrtS**4*x1) + 
     &                           (1 - 2*muu2**2)*muu4**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2))))))*HR(indj,i,tt,gen))))
     &      ))/(sqrtS**4*(1 + muu3**2 - muu4**2 - x1)*(4*muu3**2 - x1**2)*((muu2**2 + muu3**2)*sqrtS**2 + 
     &      (4*(-1 + muu1**2 - muu2**2)*muu3**4*sqrtS**6 + (-1 + muu1**2 - muu2**2)*(-1 + muu4**2)*sqrtS**6*x1**2 - 
     &         4*p1k1*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) + 
     &         (1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt(((-1 + muu3**2)**2 - 2*(1 + muu3**2)*muu4**2 + muu4**4)*sqrtS**4)*x1*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         (-1 + muu1**2 - muu2**2)*muu3**2*sqrtS**6*(-4 + 4*muu4**2 + x1**2))/(2.*sqrtS**4*(4*muu3**2 - x1**2)) - (MSf2(i,itypeU,gen)-cI*MSf(i,itypeU,gen)*WSf(i,itypeU,gen)))*
     &    ((muu1**2 + muu3**2)*sqrtS**2 - 2*(((1 + muu1**2 - muu2**2)*sqrtS**2)/4. + ((1 + muu1**2 - muu2**2)*muu3**2*sqrtS**2)/4. - ((1 + muu1**2 - muu2**2)*muu4**2*sqrtS**2)/4. + 
     &         (p1k1*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4))/Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)) - 
     &         ((1 + muu1**2 - muu2**2)*sqrtS**2*Sqrt((muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*sqrtS**4)*x1)/(4.*Sqrt(sqrtS**4*(-4*muu3**2 + x1**2)))) - (MSf2(j,itypeT,gen)+cI*MSf(j,itypeT,gen)*WSf(j,itypeT,gen))
     &      ))   )

          endif

	  enddo


	enddo

	end





