************************************************

	double precision function NeuQ2qx_M2TreeGAred()

	implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	

	integer i,j
	complex*16 DenSmZ, DenSmH(4)
	double precision colorfactor
	double precision DenSmZ1, DenSmH1(4), M2test
	double precision treeSS, treeTT, treeST, treeTS

        NeuQ2qx_M2TreeGAred	= 0d0
        treeSS = 0d0
        treeTT = 0d0
        treeST = 0d0
        treeTS = 0d0
        

! ==================================================================================================================================	
!                             ---                     gluon or photon final state           ---
! ==================================================================================================================================

	
! --- 	setting the specific couplings B,C for gluon or photon final state   ---
	   if (xtype.eq.1) then
           colorfactor = 4./3.	 
           else if (xtype.eq.2) then
           colorfactor = 1.0
           endif

!	print *, "Tree_before=", NeuQ2qx_M2TreeGAred
	
! --- Mss for gluon ---	
	treeSS = colorfactor*NC*((2*B(xtype,pstype)**2*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1red*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1red**2*(Mf(pstype,psgen)/sqrtS) +
     & m2red**2*(Mf(pstype,psgen)/sqrtS) + 2*m3red**2*(Mf(pstype,psgen)/sqrtS) - 2*m3red*(Mf(pstype,psgen)/sqrtS)**2 - 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*tred - (Mf(pstype,psgen)/sqrtS)*ured) + 
     & dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1red**4 - m2red**4 + m3red**2*(Mf(pstype,psgen)/sqrtS)**2 + m3red**2*1d0 - 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 -
     & (Mf(pstype,psgen)/sqrtS)**2*tred + m1red**2*(2*m3red**2 - 4*m3red*(Mf(pstype,psgen)/sqrtS) + (Mf(pstype,psgen)/sqrtS)**2 - tred - ured) - 1d0*ured + m2red**2*(-2*m3red**2 + 4*m3red*(Mf(pstype,psgen)/sqrtS) +
     & 1d0 + tred + ured))) + AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1red*(m1red**2*(Mf(pstype,psgen)/sqrtS) +
     & m2red**2*(Mf(pstype,psgen)/sqrtS) + 2*m3red**2*(Mf(pstype,psgen)/sqrtS) - 2*m3red*(Mf(pstype,psgen)/sqrtS)**2 - 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*tred - (Mf(pstype,psgen)/sqrtS)*ured) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1red**4 - m2red**4 + m3red**2*(Mf(pstype,psgen)/sqrtS)**2 + m3red**2*1d0 - 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 - (Mf(pstype,psgen)/sqrtS)**2*tred +
     & m1red**2*(2*m3red**2 - 4*m3red*(Mf(pstype,psgen)/sqrtS) + (Mf(pstype,psgen)/sqrtS)**2 - tred - ured) - 1d0*ured + m2red**2*(-2*m3red**2 + 4*m3red*(Mf(pstype,psgen)/sqrtS) + 1d0 + tred + ured)))))/(((Mf(pstype,psgen)/sqrtS)**2 - 1d0 + (cI*(Mf(pstype,psgen)/sqrtS)*(width(pstype)/sqrtS)))**2))
   
        NeuQ2qx_M2TreeGAred = NeuQ2qx_M2TreeGAred + treeSS
        print *, "Tree_Mss=", treeSS
!        print *, "Tree_M2Tree=", NeuQ2qx_M2TreeGAred

	
        
! --- Mtt for gluon ---
	treeTT = colorfactor*NC*(-1.*C(xtype,pstype)**2*(AL(chimeigen,ptmeigen,pttype,ptgen)*(2*m1red*m3red*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen)) + 
     & dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(m1red**2 + m3red**2 - tred)) + 
     & AR(chimeigen,ptmeigen,pttype,ptgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1red*m3red +
     & dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(m1red**2 + m3red**2 - tred)))*(2*m2red**2 - m4red**2 + 2*tred))/(((MSf(ptmeigen,pttype,ptgen)/sqrtS)**2 - tred)**2)
	
     	NeuQ2qx_M2TreeGAred = NeuQ2qx_M2TreeGAred + treeTT
   	print *, "Tree_Mtt=", treeTT
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeGAred
     	
 
! --- Mst for gluon ---
	treeST = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1red*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-3*m1red**2*m3red -
     & 2*m3red*m4red**2 + 3*m3red**2*(Mf(pstype,psgen)/sqrtS) + m4red**2*(Mf(pstype,psgen)/sqrtS) + m2red**2*(m3red + 2*(Mf(pstype,psgen)/sqrtS)) + 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*1d0 + m3red*tred + m3red*ured -
     & 2*(Mf(pstype,psgen)/sqrtS)*ured) + dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3red*m4red**2*(Mf(pstype,psgen)/sqrtS) + m3red**2*1d0 + m4red**2*1d0 +
     & 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 - 1d0**2 + m3red**2*tred + m4red**2*tred - tred**2 + m1red**2*(8*m2red**2 + 6*m3red*(m3red - (Mf(pstype,psgen)/sqrtS)) - 1d0 + tred - 5*ured) +
     & m2red**2*(2*m3red*(m3red - 2*(Mf(pstype,psgen)/sqrtS)) + 1d0 - 3*tred - ured) - (m3red**2 + m4red**2 - 2*m3red*(Mf(pstype,psgen)/sqrtS))*ured + ured**2)) +
     & AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1red*(-3*m1red**2*m3red - 2*m3red*m4red**2 +
     & 3*m3red**2*(Mf(pstype,psgen)/sqrtS) + m4red**2*(Mf(pstype,psgen)/sqrtS) + m2red**2*(m3red + 2*(Mf(pstype,psgen)/sqrtS)) + 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*1d0 + m3red*tred + m3red*ured -
     & 2*(Mf(pstype,psgen)/sqrtS)*ured) + dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3red*m4red**2*(Mf(pstype,psgen)/sqrtS) + m3red**2*1d0 + m4red**2*1d0 +
     & 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 - 1d0**2 + m3red**2*tred + m4red**2*tred - tred**2 + m1red**2*(8*m2red**2 + 6*m3red*(m3red - (Mf(pstype,psgen)/sqrtS)) - 1d0 + tred - 5*ured) +
     & m2red**2*(2*m3red*(m3red - 2*(Mf(pstype,psgen)/sqrtS)) + 1d0 - 3*tred - ured) - (m3red**2 + m4red**2 - 2*m3red*(Mf(pstype,psgen)/sqrtS))*ured + ured**2))))/(2.*((Mf(pstype,psgen)/sqrtS)**2 -
     & 1d0 + (cI*(Mf(pstype,psgen)/sqrtS)*(width(pstype)/sqrtS)))*((MSf(ptmeigen,pttype,ptgen)/sqrtS)**2 - tred))) 

        NeuQ2qx_M2TreeGAred = NeuQ2qx_M2TreeGAred + treeST
   	print *, "Tree_Mst=", treeST
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeGAred
     
       
! --- Mts for gluon ---
       treeTS = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2*m1red*AR(chimeigen,ptmeigen,pttype,ptgen)*(-3*m1red**2*m3red -
     & 2*m3red*m4red**2 + 3*m3red**2*(Mf(pstype,psgen)/sqrtS) + m4red**2*(Mf(pstype,psgen)/sqrtS) + m2red**2*(m3red + 2*(Mf(pstype,psgen)/sqrtS)) + 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*1d0 + m3red*tred + m3red*ured -
     & 2*(Mf(pstype,psgen)/sqrtS)*ured) + AL(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3red*m4red**2*(Mf(pstype,psgen)/sqrtS) + m3red**2*1d0 + m4red**2*1d0 + 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 -
     & 1d0**2 + m3red**2*tred + m4red**2*tred - tred**2 + m1red**2*(8*m2red**2 + 6*m3red*(m3red - (Mf(pstype,psgen)/sqrtS)) - 1d0 + tred - 5*ured) + m2red**2*(2*m3red*(m3red - 2*(Mf(pstype,psgen)/sqrtS)) +
     & 1d0 - 3*tred - ured) - (m3red**2 + m4red**2 - 2*m3red*(Mf(pstype,psgen)/sqrtS))*ured + ured**2)) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(chimeigen,ptmeigen,pttype,ptgen)*m1red*(-3*m1red**2*m3red - 2*m3red*m4red**2 +
     & 3*m3red**2*(Mf(pstype,psgen)/sqrtS) + m4red**2*(Mf(pstype,psgen)/sqrtS) + m2red**2*(m3red + 2*(Mf(pstype,psgen)/sqrtS)) + 2*m3red*1d0 - (Mf(pstype,psgen)/sqrtS)*1d0 + m3red*tred + m3red*ured -
     & 2*(Mf(pstype,psgen)/sqrtS)*ured) + AR(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3red*m4red**2*(Mf(pstype,psgen)/sqrtS) + m3red**2*1d0 + m4red**2*1d0 + 4*m3red*(Mf(pstype,psgen)/sqrtS)*1d0 -
     & 1d0**2 + m3red**2*tred + m4red**2*tred - tred**2 + m1red**2*(8*m2red**2 + 6*m3red*(m3red - (Mf(pstype,psgen)/sqrtS)) - 1d0 + tred - 5*ured) +  m2red**2*(2*m3red*(m3red - 2*(Mf(pstype,psgen)/sqrtS)) +
     & 1d0 - 3*tred - ured) - (m3red**2 + m4red**2 - 2*m3red*(Mf(pstype,psgen)/sqrtS))*ured + ured**2))))/(2.*((Mf(pstype,psgen)/sqrtS)**2 - 1d0 + (cI*(Mf(pstype,psgen)/sqrtS)*(width(pstype)/sqrtS)))*((MSf(ptmeigen,pttype,ptgen)/sqrtS)**2 - tred)))

        NeuQ2qx_M2TreeGAred = NeuQ2qx_M2TreeGAred + treeTS
    	print *, "Tree_Mts=", treeTS
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeGAred

	
       end
