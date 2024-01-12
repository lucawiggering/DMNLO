
***************** Fermion charges ************************

	double precision function ef(type)
	
	integer type
	
	if (type.eq.1) then
	
	 ef = 0d0
	 
	elseif (type.eq.2) then
	
	 ef = -1d0
	 
	elseif (type.eq.3) then
	
	 ef = 2/3d0
	 
	elseif (type.eq.4) then
	
	 ef = -1/3d0
	 
	else 
	endif
	end

***************** Fermion isospin ************************

	double precision function I3f(type)

	integer type
	
	if ((type.eq.1).or.(type.eq.3)) then
	
  	 I3f = 1/2d0
  	 
	 else
	 
  	 I3f = -1/2d0
  	 
	endif	
	end

	
***************** Yukawa coupling ***********************

	double precision function hf(type, gen)
#include "../DMNLO_Model.h"

	integer type, gen
	
	if ((type.eq.1).or.(type.eq.3)) then
	
  	 hf = EL*Mf(type,gen)/(sqrt2*MW*SW*SB) 
  	 
	else
	
  	 hf = EL*Mf(type,gen)/(sqrt2*MW*SW*CB) 
  	 
	endif	
	end

***********************************************

	double precision function kappa(sqs,m1,m2)

	double precision sqs, m1, m2

	 kappa = dsqrt((sqs-m1-m2)**2 - 4d0*m1*m2)
	 
	end

************************************************

	double complex function ckappa(sqs,m1,m2)

	double precision sqs, m1, m2

	 ckappa = cdsqrt(dcmplx((sqs-m1-m2)**2 - 4d0*m1*m2))

	end

************************************************

	double precision function Krondelta(i,j)

	integer i,j

	if (i.eq.j) then
	 Krondelta = 1d0 	
	else
	 Krondelta = 0d0 	
	endif

	end

************************************************
! ======= Anti-symetric 2x2 matrix ========

	double precision function Epsij(s,t)

	integer s,t
	if (s.eq.t) then
	 	Epsij =  0d0 	
	else if ((s.eq.1).and.(t.eq.2)) then
	 	Epsij =  1d0 	
	else if ((s.eq.2).and.(t.eq.1)) then
	 	Epsij = -1d0 	
	endif

	return
	end

