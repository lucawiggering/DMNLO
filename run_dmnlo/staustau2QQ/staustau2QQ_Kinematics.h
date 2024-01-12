	double precision muu1,muu2,muu3,muu4,muu5
	double precision x1,x2,tman,tred,uman,ured,xx1,xx2,xx3,jacobi
	double precision flux,sqrtS, DenSV(2)
	double precision p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3,p1p2

	integer isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr
	integer ftt1, ftt2, fgen1, fgen2
	integer imin, imax, jmin, jmax, bmin, bmax, cmin, cmax, pmin, pmax, qmin, qmax

	common /mass/ muu1,muu2,muu3,muu4,muu5
	common /kinematics/ flux,sqrtS, DenSV
	common /var23/ x1, x2,tman,tred,uman,ured,xx1,xx2,xx3,jacobi
	common /var2to3/ p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3,p1p2


	common /typegg/ isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr
	common /typegg/ ftt1, ftt2, fgen1, fgen2
	common /typegg/ imin, imax, jmin, jmax, bmin, bmax, cmin, cmax, pmin, pmax, qmin, qmax
