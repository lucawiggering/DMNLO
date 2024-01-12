************************************************

	double precision function ChiChi2QQ_M2TreeM0()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j
	complex*16 DenSmZ, DenSmH(4)
	double precision DenSmZ1, DenSmH1(4)

	ChiChi2QQ_M2TreeM0 = 0d0
	
c	wZ = 0d0
c	wHk0(1) = 0d0
c	wHk0(2) = 0d0
c	wHk0(3) = 0d0
c	wHk0(4) = 0d0

************ Zboson contr.

	DenSmZ = (sqrtS**2 - MZ2) + cI*MZ*wZ
	DenSmZ1 = DenSmZ*dconjg(DenSmZ)/(sqrtS**2 - MZ2)
	

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (16*AL**2*NC*(BL*BR*(1 - 6*muu1**2)*muu3**2*sqrtS**4 + ((BL**2 + BR**2)*
     -         ((1 - 2*muu3**2 + 2*(muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2)))*sqrtS**4 - 
     -          2*(-1 + 2*muu1**2 + 2*muu3**2)*sqrtS**2*tman + 2*tman**2))/4d0))/(DenSmZ*dconjg(DenSmZ))


************ Higgs contr.

	do i = 1, 4
	  DenSmH(i) = (sqrtS**2 - MHk0(i)**2) + cI*MHk0(i)*wHk0(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk0(i)**2)
	enddo

	do i = 1, 4

c	  write(*,*) 'Couplings = ',DL(i), DR(i), DLM(i), DRM(i)
	
	  do j = 1, 4

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (NC*sqrtS**4*(-(CL(i)*CL(j)) - CR(i)*CR(j) + 2*muu1**2*(CL(i) + CR(i))*(CL(j) + CR(j)))*
     -    (-(DLM(i)*DLM(j)) - DRM(i)*DRM(j) + 2*muu3**2*(DLM(i) + DRM(i))*(DLM(j) + DRM(j))))/
     -  (DenSmH(i)*dconjg(DenSmH(j)))

	  enddo

************ ZHiggs contr.

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (4*AL*(BL - BR)*muu1*muu3*NC*sqrtS**4*(CL(i) - CR(i))*(DLM(i) - DRM(i)))/
     -                           (dreal(DenSmZ*dconjg(DenSmH(i))))

	enddo

c	print*,'ChiChi2QQ_M2TreeM0', ChiChi2QQ_M2TreeM0

	do i = 1, 2
	  do j = 1, 2

************ T-channel contr.
	

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (4*NC*(muu1*muu3*sqrtS**2*(ELcM(j,tt,gen)*ERM(i,tt,gen) + ELcM(i,tt,gen)*ERM(j,tt,gen)) + 
     -      ((-((-1 + muu1**2 + muu3**2)*sqrtS**2) + tman)*(ELcM(i,tt,gen)*ELcM(j,tt,gen) + ERM(i,tt,gen)*ERM(j,tt,gen)))/2.)*
     -    (muu1*muu3*sqrtS**2*(FLM(j,tt,gen)*FRM(i,tt,gen) + FLM(i,tt,gen)*FRM(j,tt,gen)) + 
     -      ((-((-1 + muu1**2 + muu3**2)*sqrtS**2) + tman)*(FLM(i,tt,gen)*FLM(j,tt,gen) + FRM(i,tt,gen)*FRM(j,tt,gen)))/2.))/
     -  ((sqrtS**2 - 2*muu1**2*sqrtS**2 - 2*muu3**2*sqrtS**2 + tman + MSf(i,tt,gen)**2)*
     -    (sqrtS**2 - 2*muu1**2*sqrtS**2 - 2*muu3**2*sqrtS**2 + tman + MSf(j,tt,gen)**2))

	  enddo
	enddo

	do i = 1, 2
	  do j = 1, 2

************ U-channel contr.

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (NC*(-2*muu1*muu3*sqrtS**2*(GLM(j,tt,gen)*GRM(i,tt,gen) + GLM(i,tt,gen)*GRM(j,tt,gen)) - 
     -      muu3**2*sqrtS**2*(GLM(i,tt,gen)*GLM(j,tt,gen) + GRM(i,tt,gen)*GRM(j,tt,gen)) + 
     -      (-(muu1**2*sqrtS**2) + tman)*(GLM(i,tt,gen)*GLM(j,tt,gen) + GRM(i,tt,gen)*GRM(j,tt,gen)))*
     -    (-2*muu1*muu3*sqrtS**2*(HLM(j,tt,gen)*HRM(i,tt,gen) + HLM(i,tt,gen)*HRM(j,tt,gen)) - 
     -      muu3**2*sqrtS**2*(HLM(i,tt,gen)*HLM(j,tt,gen) + HRM(i,tt,gen)*HRM(j,tt,gen)) + 
     -      (-(muu1**2*sqrtS**2) + tman)*(HLM(i,tt,gen)*HLM(j,tt,gen) + HRM(i,tt,gen)*HRM(j,tt,gen))))/
     -                  ((-tman + MSf(i,tt,gen)**2)*(-tman + MSf(j,tt,gen)**2))

	  enddo
	enddo

	do i = 1, 2
	
	  do j = 1, 2

************ TU-channel contr. (minus for a wrong sign in FORM)

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 - (2*NC*(muu1**4*sqrtS**4*(ELcM(i,tt,gen)*FLM(i,tt,gen)*GLM(j,tt,gen)*HLM(j,tt,gen) + 
     -       ERM(i,tt,gen)*FRM(i,tt,gen)*GRM(j,tt,gen)*HRM(j,tt,gen)) + 
     -      muu3**4*sqrtS**4*(ELcM(i,tt,gen)*FLM(i,tt,gen)*GLM(j,tt,gen)*HLM(j,tt,gen) + ERM(i,tt,gen)*FRM(i,tt,gen)*GRM(j,tt,gen)*HRM(j,tt,gen)) + 
     -      tman*(sqrtS**2 + tman)*(ELcM(i,tt,gen)*FLM(i,tt,gen)*GLM(j,tt,gen)*HLM(j,tt,gen) + ERM(i,tt,gen)*FRM(i,tt,gen)*GRM(j,tt,gen)*HRM(j,tt,gen)) + 
     -      muu1**2*sqrtS**2*(ELcM(i,tt,gen)*(-2*tman*FLM(i,tt,gen)*GLM(j,tt,gen) + sqrtS**2*FRM(i,tt,gen)*GRM(j,tt,gen))*HLM(j,tt,gen) + 
     -         ERM(i,tt,gen)*(sqrtS**2*FLM(i,tt,gen)*GLM(j,tt,gen) - 2*tman*FRM(i,tt,gen)*GRM(j,tt,gen))*HRM(j,tt,gen)) + 
     -      muu1*muu3**3*sqrtS**4*(ERM(i,tt,gen)*(FRM(i,tt,gen)*(GRM(j,tt,gen)*HLM(j,tt,gen) + GLM(j,tt,gen)*HRM(j,tt,gen)) - 
     -            FLM(i,tt,gen)*(GLM(j,tt,gen)*HLM(j,tt,gen) + GRM(j,tt,gen)*HRM(j,tt,gen))) + 
     -         ELcM(i,tt,gen)*(FLM(i,tt,gen)*(GRM(j,tt,gen)*HLM(j,tt,gen) + GLM(j,tt,gen)*HRM(j,tt,gen)) - 
     -            FRM(i,tt,gen)*(GLM(j,tt,gen)*HLM(j,tt,gen) + GRM(j,tt,gen)*HRM(j,tt,gen)))) - 
     -      muu1*muu3*sqrtS**2*(ERM(i,tt,gen)*((-(muu1**2*sqrtS**2) + tman)*FRM(i,tt,gen)*(GRM(j,tt,gen)*HLM(j,tt,gen) + GLM(j,tt,gen)*HRM(j,tt,gen)) + 
     -            (-sqrtS**2 + muu1**2*sqrtS**2 - tman)*FLM(i,tt,gen)*(GLM(j,tt,gen)*HLM(j,tt,gen) + GRM(j,tt,gen)*HRM(j,tt,gen))) + 
     -         ELcM(i,tt,gen)*((-(muu1**2*sqrtS**2) + tman)*FLM(i,tt,gen)*(GRM(j,tt,gen)*HLM(j,tt,gen) + GLM(j,tt,gen)*HRM(j,tt,gen)) + 
     -            (-sqrtS**2 + muu1**2*sqrtS**2 - tman)*FRM(i,tt,gen)*(GLM(j,tt,gen)*HLM(j,tt,gen) + GRM(j,tt,gen)*HRM(j,tt,gen)))) + 
     -      muu3**2*sqrtS**2*(-2*tman*ELcM(i,tt,gen)*FLM(i,tt,gen)*GLM(j,tt,gen)*HLM(j,tt,gen) + 
     -                    sqrtS**2*ERM(i,tt,gen)*FLM(i,tt,gen)*GRM(j,tt,gen)*HLM(j,tt,gen) + 
     -         sqrtS**2*ELcM(i,tt,gen)*FRM(i,tt,gen)*GLM(j,tt,gen)*HRM(j,tt,gen) - 2*tman*ERM(i,tt,gen)*FRM(i,tt,gen)*GRM(j,tt,gen)*HRM(j,tt,gen) - 
     -         2*muu1**2*sqrtS**2*(ERM(i,tt,gen)*(-(FRM(i,tt,gen)*GLM(j,tt,gen)*HLM(j,tt,gen)) + FLM(i,tt,gen)*GRM(j,tt,gen)*HLM(j,tt,gen) + 
     -               FLM(i,tt,gen)*GLM(j,tt,gen)*HRM(j,tt,gen)) + ELcM(i,tt,gen)*
     -             (FRM(i,tt,gen)*GRM(j,tt,gen)*HLM(j,tt,gen) + FRM(i,tt,gen)*GLM(j,tt,gen)*HRM(j,tt,gen) - FLM(i,tt,gen)*GRM(j,tt,gen)*HRM(j,tt,gen))))))/
     -  ((MSf(j,tt,gen)**2 - tman)*(MSf(i,tt,gen)**2 + sqrtS**2 - 2*muu1**2*sqrtS**2 - 2*muu3**2*sqrtS**2 + tman))

	  enddo
	enddo

	do i = 1, 2

************ TS-channel Zboson contr.

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (4*AL*NC*(-(muu3**4*sqrtS**4*(BL*ERM(i,tt,gen)*FLM(i,tt,gen) - BR*ELcM(i,tt,gen)*FRM(i,tt,gen))) - 
     -      (muu1**4*sqrtS**4 + (sqrtS**2 + tman)**2 - muu1**2*sqrtS**2*(3*sqrtS**2 + 2*tman))*
     -        (BL*ERM(i,tt,gen)*FLM(i,tt,gen) - BR*ELcM(i,tt,gen)*FRM(i,tt,gen)) + 
     -      3*(BL - BR)*muu1*muu3**3*sqrtS**4*(ELcM(i,tt,gen)*FLM(i,tt,gen) + ERM(i,tt,gen)*FRM(i,tt,gen)) + 
     -      (BL - BR)*muu1*muu3*sqrtS**2*((-2 + 3*muu1**2)*sqrtS**2 - 3*tman)*(ELcM(i,tt,gen)*FLM(i,tt,gen) + ERM(i,tt,gen)*FRM(i,tt,gen)) - 
     -      muu3**2*sqrtS**2*(BL*(4*muu1**2*sqrtS**2*ERM(i,tt,gen)*FLM(i,tt,gen) - 2*(sqrtS**2 + tman)*ERM(i,tt,gen)*FLM(i,tt,gen) + 
     -            (-1 + 6*muu1**2)*sqrtS**2*ELcM(i,tt,gen)*FRM(i,tt,gen)) + 
     -         BR*((1 - 6*muu1**2)*sqrtS**2*ERM(i,tt,gen)*FLM(i,tt,gen) + 2*(sqrtS**2 - 2*muu1**2*sqrtS**2 + tman)*ELcM(i,tt,gen)*FRM(i,tt,gen)))))/
     -  (-DenSmZ1*(sqrtS**2 - 2*muu1**2*sqrtS**2 - 2*muu3**2*sqrtS**2 + tman + MSf(i,tt,gen)**2))

	  do j = 1, 4

************ TS-channel Hboson contr.

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 + (-2*NC*sqrtS**2*(-(muu1*muu3**3*sqrtS**2*(CL(j) + CR(j))*(DLM(j) + DRM(j))*
     -                (ERM(i,tt,gen)*FLM(i,tt,gen) + ELcM(i,tt,gen)*FRM(i,tt,gen))) + 
     -      muu1*muu3*(sqrtS**2*CL(j)*DLM(j) + tman*CL(j)*DLM(j) + tman*CR(j)*DLM(j) + tman*CL(j)*DRM(j) + sqrtS**2*CR(j)*DRM(j) + tman*CR(j)*DRM(j) - 
     -         muu1**2*sqrtS**2*(CL(j) + CR(j))*(DLM(j) + DRM(j)))*(ERM(i,tt,gen)*FLM(i,tt,gen) + ELcM(i,tt,gen)*FRM(i,tt,gen)) + 
     -      (sqrtS**2 + tman)*(CL(j)*DLM(j)*ELcM(i,tt,gen)*FLM(i,tt,gen) + CR(j)*DRM(j)*ERM(i,tt,gen)*FRM(i,tt,gen)) - 
     -      muu1**2*sqrtS**2*((2*CL(j) + CR(j))*DLM(j)*ELcM(i,tt,gen)*FLM(i,tt,gen) + (CL(j) + 2*CR(j))*DRM(j)*ERM(i,tt,gen)*FRM(i,tt,gen)) + 
     -      muu3**2*sqrtS**2*(-2*CL(j)*DLM(j)*ELcM(i,tt,gen)*FLM(i,tt,gen) - CL(j)*DRM(j)*ELcM(i,tt,gen)*FLM(i,tt,gen) - 
     -       CR(j)*DLM(j)*ERM(i,tt,gen)*FRM(i,tt,gen) - 2*CR(j)*DRM(j)*ERM(i,tt,gen)*FRM(i,tt,gen) + 2*muu1**2*(CL(j) + 
     -       CR(j))*(DLM(j) + DRM(j))*(ELcM(i,tt,gen)*FLM(i,tt,gen) + ERM(i,tt,gen)*FRM(i,tt,gen)))))/
     -  (-DenSmH1(j)*(sqrtS**2 - 2*muu1**2*sqrtS**2 - 2*muu3**2*sqrtS**2 + tman + MSf(i,tt,gen)**2))

	  enddo
	enddo

	do i = 1, 2

************ US-channel Zboson contr. (minus for a wrong sign in FORM)

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 - (4*AL*NC*(-(muu3**4*sqrtS**4*(BR*GRM(i,tt,gen)*HLM(i,tt,gen) - BL*GLM(i,tt,gen)*HRM(i,tt,gen))) - 
     -      (-(muu1**2*sqrtS**4) + (-(muu1**2*sqrtS**2) + tman)**2)*(BR*GRM(i,tt,gen)*HLM(i,tt,gen) - BL*GLM(i,tt,gen)*HRM(i,tt,gen)) + 
     -      3*(BL - BR)*muu1*muu3**3*sqrtS**4*(GLM(i,tt,gen)*HLM(i,tt,gen) + GRM(i,tt,gen)*HRM(i,tt,gen)) + 
     -      (-BL + BR)*muu1*muu3*sqrtS**2*(sqrtS**2 - 3*muu1**2*sqrtS**2 + 3*tman)*(GLM(i,tt,gen)*HLM(i,tt,gen) + GRM(i,tt,gen)*HRM(i,tt,gen)) + 
     -      muu3**2*sqrtS**2*(-(BL*sqrtS**2*GRM(i,tt,gen)*HLM(i,tt,gen)) + 2*BR*tman*GRM(i,tt,gen)*HLM(i,tt,gen) + BR*sqrtS**2*GLM(i,tt,gen)*HRM(i,tt,gen) - 
     -         2*BL*tman*GLM(i,tt,gen)*HRM(i,tt,gen) + 2*muu1**2*sqrtS**2*
     -          (3*BL*GRM(i,tt,gen)*HLM(i,tt,gen) - 2*BR*GRM(i,tt,gen)*HLM(i,tt,gen) + 2*BL*GLM(i,tt,gen)*HRM(i,tt,gen) - 3*BR*GLM(i,tt,gen)*HRM(i,tt,gen)))))/
     -  (-DenSmZ1*(-tman + MSf(i,tt,gen)**2))

	  do j = 1, 4

************ US-channel Hboson contr. (minus for a wrong sign in FORM)

	ChiChi2QQ_M2TreeM0 = ChiChi2QQ_M2TreeM0 - (2*NC*sqrtS**2*(-(tman*CL(j)*DLM(j)*GLM(i,tt,gen)*HLM(i,tt,gen)) - muu1**2*sqrtS**2*CR(j)*DLM(j)*GLM(i,tt,gen)*HLM(i,tt,gen) - 
     -      muu1**2*sqrtS**2*CL(j)*DRM(j)*GRM(i,tt,gen)*HRM(i,tt,gen) - tman*CR(j)*DRM(j)*GRM(i,tt,gen)*HRM(i,tt,gen) + 
     -      muu1*muu3**3*sqrtS**2*(CL(j) + CR(j))*(DLM(j) + DRM(j))*(GRM(i,tt,gen)*HLM(i,tt,gen) + GLM(i,tt,gen)*HRM(i,tt,gen)) + 
     -      muu1*muu3*(muu1**2*sqrtS**2*(CL(j) + CR(j))*(DLM(j) + DRM(j)) - tman*(CL(j) + CR(j))*(DLM(j) + DRM(j)) - sqrtS**2*(CR(j)*DLM(j) + CL(j)*DRM(j)))*
     -       (GRM(i,tt,gen)*HLM(i,tt,gen) + GLM(i,tt,gen)*HRM(i,tt,gen)) + 
     -      muu3**2*sqrtS**2*(-(CL(j)*DRM(j)*GLM(i,tt,gen)*HLM(i,tt,gen)) - CR(j)*DLM(j)*GRM(i,tt,gen)*HRM(i,tt,gen) + 
     -         2*muu1**2*(CL(j) + CR(j))*(DLM(j) + DRM(j))*(GLM(i,tt,gen)*HLM(i,tt,gen) + GRM(i,tt,gen)*HRM(i,tt,gen)))))/
     -  (-DenSmH1(j)*(-tman + MSf(i,tt,gen)**2))


	  enddo
	enddo



	end
