      double precision muu1,muu2,muu3,muu4,muu5,MGlr
      double precision jacobi,xx1,xx2,xx3
      double precision kakb,k1k2,k1k3,k2k3,kak1,kak2,kak3,kbk1,kbk2,kbk3
      double precision pap1, pap2, pap3, pbp1, pbp2, pbp3, p1p2, p1p3, p2p3, papb
      double precision x1,x2,p1k1,p1k2,tman,tred,uman,ured
      double precision flux,sqrtS,sman
      double precision x0,kappaDip


      integer fluxflag

      common /kinematics/ flux,sqrtS,sman
      common /reducedmass2/ muu1, muu2, muu3, muu4, muu5, MGlr 
      common /genvar23/ jacobi,xx1,xx2,xx3
      common /genvar23/ kakb,k1k2,k1k3,k2k3,kak1,kak2,kak3,kbk1,kbk2,kbk3      
      common /var23/ pap1, pap2, pap3, pbp1, pbp2, pbp3, p1p2, p1p3, p2p3, papb
      common /var23/ x1, x2, p1k1, p1k2,tman,tred,uman,ured
      common /typegg/ fluxflag
      common /dipole/ x0,kappaDip
      