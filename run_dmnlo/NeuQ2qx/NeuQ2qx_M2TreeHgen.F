************************************************

	double precision function NeuQ2qx_M2TreeHgen()

	implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	external myCget, myBget

	
	complex*16 DenSmZ, DenSmH(4)
	double precision colorfactor
	double precision DenSmZ1, DenSmH1(4), M2test
	double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
	
	double precision MssH, MttH, MuuH, MstH, MtsH, MsuH, MusH, MutH, MtuH
	
	double complex A1vL,A1vR,A1sL,A1sR,A2sL,A2sR,A3vL,A3vR,A3sL,A3sR
	double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc
	double complex FvL,FvR,FsL,FsR
	double complex FLc,FRc
	double complex Gscalar,Gc
	double complex HtreeL,HtreeR,HLc,HRc
	
	double complex PvL,PvR,PsL,PsR,P

	double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
	
	integer sumU, sumUc, sumT, sumTc,j
	integer sumUend, sumUcend, sumTend, sumTcend
	double complex Gl(2)

        NeuQ2qx_M2TreeHgen	= 0d0
        treeSS = 0d0
        treeTT = 0d0
        treeUU = 0d0
        treeST = 0d0
        treeTS = 0d0
        treeSU = 0d0
        treeUS = 0d0
        treeUT = 0d0
        treeTU = 0d0
        
        sumU = 0
        sumUc = 0
        sumT = 0
        sumTc = 0

        
            
        if ((xtype.eq.5).or.(xtype.eq.6)) then
        
           sumTend = 1
           sumTcend = 1
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.7) then
        
           sumTend = 0
           sumTcend = 0
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.8) then
               
           sumTend = 1
           sumTcend = 1
           sumUend = 1
           sumUcend = 1
           
         endif
                

! ==================================================================================================================================	
!                             ---                     Higgs final states           ---
! ==================================================================================================================================


	
! --- Mss  ---	

	
	A1vL = 0d0
	A1vR = 0d0
	A1sL = AL(chimeigen,sfmeigen, sftype,sfgen)
	A1sR = AR(chimeigen,sfmeigen, sftype,sfgen)
	PvL = 1d0
	PvR = 1d0
	PsL = Mf(pstype,psgen)
	PsR = Mf(pstype,psgen)
	FvL = 0d0
	FvR = 0d0
	FsL = FL(xtype,pstype,psgen)
	FsR = FR(xtype,pstype,psgen)
	
	A1Lc = dconjg(AL(chimeigen,sfmeigen, sftype,sfgen))
	A1Rc = dconjg(AR(chimeigen,sfmeigen, sftype,sfgen))
	FLc = dconjg(FL(xtype,pstype,psgen))
	FRc = dconjg(FR(xtype,pstype,psgen))
	
	mq1 = Mf(pstype,psgen) 
	mq1c = Mf(pstype,psgen)
	wquark = 0d0
	
	treeSS = treeSS + NC*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)
	
     
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeSS


	
        
        
! --- Mtt  ---

	do sumT = 0,sumTend
	   do sumTc = 0,sumTcend
	   
	      A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	      A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
	      P = 1d0
	      Gscalar = G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)
	   
	      A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	      A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	      Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))
	   
	      mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
	      mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
	      
	   
	      treeTT = treeTT + NC*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)

           enddo
        enddo
        
     	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeTT

     	
     	
     	
     	
! --- Muu  ---
	do sumU = 0,sumUend
	   do sumUc = 0,sumUcend
	   
	   A3vL = 0d0
	   A3vR = 0d0
	   A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
	   A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
	   HtreeL = HL(chimeigen,pumeigen+sumU,xtype)
	   HtreeR = HR(chimeigen,pumeigen+sumU,xtype)
	   A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	   A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	   HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	   HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))
	   mq3 =  MNeuCha(pumeigen+sumU)
	   mq3c = MNeuCha(pumeigen+sumUc)
	   
	treeUU = treeUU + NC*MuuH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A3Lc,A3Rc,HLc,HRc,mq3,mq3c)
	
           enddo
        enddo
        
     	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeUU

     	     	
! --- Mst  ---

	do sumTc = 0,sumTcend
	
	A1vL = 0d0
	A1vR = 0d0
	A1sL = AL(chimeigen,sfmeigen, sftype,sfgen)
	A1sR = AR(chimeigen,sfmeigen, sftype,sfgen)
	PvL = 1d0
	PvR = 1d0
	PsL = Mf(pstype,psgen)
	PsR = Mf(pstype,psgen)
	FvL = 0d0
	FvR = 0d0
	FsL = FL(xtype,pstype,psgen)
	FsR = FR(xtype,pstype,psgen)
	A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))
	A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen))
	Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))
	
	mq1 = Mf(pstype,psgen)
	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
	wquark = 0d0

	treeST = treeST + NC*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)
	
	enddo

        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeST

     
       
     	
     	
! --- Mts  ---
  
	do sumT = 0,sumTend
	
	A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
	P = 1d0
	Gscalar = G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)
	
	A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	FLc = dconjg(FL(xtype,pstype,psgen))
	FRc = dconjg(FR(xtype,pstype,psgen))
	
	mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
	mq1c = Mf(pstype,psgen)
	wquark = 0d0
	
	
	
	treeTS = treeTS + NC*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)

	enddo
	
	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeTS

     	
   
     	
     	
! --- Msu  ---
  
	do sumUc = 0,sumUcend
	
	A1vL = 0d0
	A1vR = 0d0
	A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
	A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)
	PvL = 1d0
	PvR = 1d0
	PsL = Mf(pstype,psgen)
	PsR = Mf(pstype,psgen)
	FvL = 0d0
	FvR = 0d0
	FsL = FL(xtype,pstype,psgen)
	FsR = FR(xtype,pstype,psgen)
	A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))
	
	mq1 = Mf(pstype,psgen)
	mq3c = MNeuCha(pumeigen+sumUc)
	wquark = 0d0
	
	treeSU = treeSU + NC*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)

	enddo
	
	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeSU

     	
     	
     	
     	
! --- Mus  ---
  
	do sumU = 0,sumUend
	
	A3vL = 0d0
	A3vR = 0d0
	A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
	A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
	HtreeL = HL(chimeigen,pumeigen+sumU,xtype)
	HtreeR = HR(chimeigen,pumeigen+sumU,xtype)
	A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	FLc = dconjg(FL(xtype,pstype,psgen))
	FRc = dconjg(FR(xtype,pstype,psgen))
	
	mq3 = MNeuCha(pumeigen+sumU)
	mq1c = Mf(pstype,psgen)
	wquark = 0d0
	
	
	treeUS = treeUS + NC*MusH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A1Lc,A1Rc,FLc,FRc,mq3,mq1c,wquark)
	
	enddo
	
	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeUS


     	
     	
     	
! --- Mtu  ---
  
	do sumT = 0,sumTend
	   do sumUc = 0,sumUcend

	   A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	   A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
	   P = 1d0
	   Gscalar = G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)
	   A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	   A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	   HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	   HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))
	   
	   mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
	   mq3c = MNeuCha(pumeigen+sumUc)
	
	
	  treeTU = treeTU + NC*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)




    	 enddo
	   enddo

	
	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeTU

     	     	
     	
     	
     	
! --- Mut  ---
  
	do sumTc = 0,sumTcend
	   do sumU = 0,sumUend
	   
	   A3vL = 0d0
	   A3vR = 0d0
	   A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
	   A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
	   HtreeL = HL(chimeigen,pumeigen+sumU,xtype)
	   HtreeR = HR(chimeigen,pumeigen+sumU,xtype)
	   A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	   A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	   Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))
	   mq3 = MNeuCha(pumeigen+sumU)
	   mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
	
	  treeUT = treeUT + NC*MutH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A2Lc,A2Rc,Gc,mq3,mq2c)
	
	   enddo
	enddo
	
	
        NeuQ2qx_M2TreeHgen = NeuQ2qx_M2TreeHgen + treeUT
	
	treeperc(1) = treeSS
	treeperc(2) = treeTT
	treeperc(3) = treeUU
	treeperc(4) = treeST + treeTS
	treeperc(5) = treeSU + treeUS
	treeperc(6) = treeTU + treeUT
	treeperc(7) = treeperc(1) + treeperc(2) + treeperc(3) + treeperc(4) + treeperc(5) + treeperc(6)
	
	
       end
