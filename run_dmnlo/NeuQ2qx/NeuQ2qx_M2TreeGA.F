************************************************

	double precision function NeuQ2qx_M2TreeGA()

	implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	

	integer i,j
	complex*16 DenSmZ, DenSmH(4)
	double precision colorfactor
	double precision DenSmZ1, DenSmH1(4), M2test
	double precision treeSS, treeTT, treeST, treeTS

        NeuQ2qx_M2TreeGA	= 0d0
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
	
! --- Mss for gluon ---	
	treeSS = colorfactor*NC*((2*B(xtype,pstype)**2*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2*Mf(pstype,psgen) +
     & m2**2*Mf(pstype,psgen) + 2*m3**2*Mf(pstype,psgen) - 2*m3*Mf(pstype,psgen)**2 - 2*m3*sman - Mf(pstype,psgen)*tman - Mf(pstype,psgen)*uman) + 
     & dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + m3**2*Mf(pstype,psgen)**2 + m3**2*sman - 4*m3*Mf(pstype,psgen)*sman -
     & Mf(pstype,psgen)**2*tman + m1**2*(2*m3**2 - 4*m3*Mf(pstype,psgen) + Mf(pstype,psgen)**2 - tman - uman) - sman*uman + m2**2*(-2*m3**2 + 4*m3*Mf(pstype,psgen) +
     & sman + tman + uman))) + AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(m1**2*Mf(pstype,psgen) +
     & m2**2*Mf(pstype,psgen) + 2*m3**2*Mf(pstype,psgen) - 2*m3*Mf(pstype,psgen)**2 - 2*m3*sman - Mf(pstype,psgen)*tman - Mf(pstype,psgen)*uman) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + m3**2*Mf(pstype,psgen)**2 + m3**2*sman - 4*m3*Mf(pstype,psgen)*sman - Mf(pstype,psgen)**2*tman +
     & m1**2*(2*m3**2 - 4*m3*Mf(pstype,psgen) + Mf(pstype,psgen)**2 - tman - uman) - sman*uman + m2**2*(-2*m3**2 + 4*m3*Mf(pstype,psgen) + sman + tman + uman)))))/((Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))**2))
   
        NeuQ2qx_M2TreeGA = NeuQ2qx_M2TreeGA + treeSS

	
        
! --- Mtt for gluon ---
	treeTT = colorfactor*NC*(-1.*C(xtype,pstype)**2*(AL(chimeigen,ptmeigen,pttype,ptgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen)) + 
     & dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(m1**2 + m3**2 - tman)) + 
     & AR(chimeigen,ptmeigen,pttype,ptgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1*m3 +
     & dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(m1**2 + m3**2 - tman)))*(2*m2**2 - m4**2 + 2*tman))/((MSf(ptmeigen,pttype,ptgen)**2 - tman)**2)
	
     	NeuQ2qx_M2TreeGA = NeuQ2qx_M2TreeGA + treeTT

! --- Mst for gluon ---
	treeST = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-3*m1**2*m3 -
     & 2*m3*m4**2 + 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sman - Mf(pstype,psgen)*sman + m3*tman + m3*uman -
     & 2*Mf(pstype,psgen)*uman) + dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sman + m4**2*sman +
     & 4*m3*Mf(pstype,psgen)*sman - sman**2 + m3**2*tman + m4**2*tman - tman**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sman + tman - 5*uman) +
     & m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) + sman - 3*tman - uman) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*uman + uman**2)) +
     & AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1*(-3*m1**2*m3 - 2*m3*m4**2 +
     & 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sman - Mf(pstype,psgen)*sman + m3*tman + m3*uman -
     & 2*Mf(pstype,psgen)*uman) + dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sman + m4**2*sman +
     & 4*m3*Mf(pstype,psgen)*sman - sman**2 + m3**2*tman + m4**2*tman - tman**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sman + tman - 5*uman) +
     & m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) + sman - 3*tman - uman) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*uman + uman**2))))/(2.*(Mf(pstype,psgen)**2 -
     & sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen,pttype,ptgen)**2 - tman))) 

        NeuQ2qx_M2TreeGA = NeuQ2qx_M2TreeGA + treeST
     
       
! --- Mts for gluon ---
       treeTS = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2*m1*AR(chimeigen,ptmeigen,pttype,ptgen)*(-3*m1**2*m3 -
     & 2*m3*m4**2 + 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sman - Mf(pstype,psgen)*sman + m3*tman + m3*uman -
     & 2*Mf(pstype,psgen)*uman) + AL(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sman + m4**2*sman + 4*m3*Mf(pstype,psgen)*sman -
     & sman**2 + m3**2*tman + m4**2*tman - tman**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sman + tman - 5*uman) + m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) +
     & sman - 3*tman - uman) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*uman + uman**2)) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(chimeigen,ptmeigen,pttype,ptgen)*m1*(-3*m1**2*m3 - 2*m3*m4**2 +
     & 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sman - Mf(pstype,psgen)*sman + m3*tman + m3*uman -
     & 2*Mf(pstype,psgen)*uman) + AR(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sman + m4**2*sman + 4*m3*Mf(pstype,psgen)*sman -
     & sman**2 + m3**2*tman + m4**2*tman - tman**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sman + tman - 5*uman) +  m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) +
     & sman - 3*tman - uman) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*uman + uman**2))))/(2.*(Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen,pttype,ptgen)**2 - tman)))

        NeuQ2qx_M2TreeGA = NeuQ2qx_M2TreeGA + treeTS

	
       end
