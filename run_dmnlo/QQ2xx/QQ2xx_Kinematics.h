	double precision muu1,muu2,muu3,muu4,muu5,MGlr
	double precision nu1, nu2, nu3, nu4
	double precision x1,x2,x3,p1k1,p1k2,tman,tred,uman,ured,cutoff,xx1,xx2,xx3,jacobi
	double precision flux,sqrtS,WeightS(6),WeightV(4)

	integer icha,indi,indj,fluxflag,total_chargeout,total_charge22,total_charge23,abort22,abort23

	common /kinematics/ flux,sqrtS,WeightS,WeightV
	common /reducedmass2/ muu1, muu2, muu3, muu4, muu5, MGlr 
	common /reducedmass3/ nu1, nu2, nu3, nu4
	common /var23/ x1, x2, x3, p1k1, p1k2,tman,tred,uman,ured,cutoff,xx1,xx2,xx3,jacobi
	common /typegg/ icha, indi,indj,fluxflag,total_chargeout,total_charge22,total_charge23,abort22,abort23
