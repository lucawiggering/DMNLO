      double precision muu1,muu2,muu3,muu4,muu5,m1
      double precision x1,x2,tman,tred,uman,ured,sman,xx1,xx2,xx3,jacobi
      double precision flux,sqrtS
      double precision x0FEISQ,x0IEISQ,x0FEISg,x0IEISg
      double precision p1p2,p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3
      double precision p1Barp2,p1Bark1,p1Bark3,p1Bark2,k1Bark2,k1Bark3,k2Bark3,p2Bark1,p2Bark3,p2Bark2,k1Bark1,k2Bark2,k3Bark3,p1Barp1,p2Barp2
      integer isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr

      common /reducedmass2/ muu1, muu2, muu3, muu4,muu5,m1
      common /var23/ x1,x2,tman,tred,uman,ured,sman,xx1,xx2,xx3,jacobi
	common /kinematics/ flux,sqrtS
      common /dipole/ x0FEISQ,x0IEISQ,x0FEISg,x0IEISg
      common /var2to3/ p1p2,p1k1,p1k2,p1k3,p2k1,p2k2,p2k3,k1k2,k1k3,k2k3
      common /forPhysicalGauge/ p1Barp2,p1Bark1,p1Bark3,p1Bark2,k1Bark2,k1Bark3,k2Bark3,p2Bark1,p2Bark3,p2Bark2,k1Bark1,k2Bark2,k3Bark3,p1Barp1,p2Barp2
      common /typegg/ isf1, itt1, igen1, isf2, itt2, igen2, fluxflag, itree, imicr
      

