	double precision muu1,muu2,muu3,muu4,muu5,m1,m2,m3,m4
	double precision muuqt(6),muupt(6),muuqu(6),muupu(6),muugl
	double precision MGau(6)
	double precision x1,x2,tman,tred,uman,ured,xx1,xx2,xx3,jacobi
	double precision p1p2,p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3
	double precision flux,sqrtS
	double precision iQCD(3),iQCDnlo(3)
	double complex vrel
	double precision x0FEIS,x0IEIS

	complex*16  DenGauT(6),DenGauU(6), DenGlT, DenGlU

	
	integer isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr
	integer ftt1, ftt2, fgen1, fgen2
	integer imin, imax, jmin, jmax, bmin, bmax, cmin, cmax, pmin, pmax, qmin, qmax
	
	common /kinematics/ flux,sqrtS
	common /reducedmass2/ muu1, muu2, muu3, muu4,muu5,m1,m2,m3,m4
	common /propmassesred/ muuqt,muupt,muuqu,muupu,muugl
	common /Gauginomass/ MGau
	common /var23/ x1, x2, tman,tred,uman,ured,xx1,xx2,xx3,jacobi
	common /var2to3/ p1p2,p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3
	common /QCDflag/ iQCD, iQCDnlo
	common /exp/ vrel
	common /dipole/ x0FEIS,x0IEIS

	common /denprop/  DenGauT,DenGauU, DenGlT, DenGlU

	common /typegg/ isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr
	common /typegg/ ftt1, ftt2, fgen1, fgen2
	common /typegg/ imin, imax, jmin, jmax, bmin, bmax, cmin, cmax, pmin, pmax, qmin, qmax
