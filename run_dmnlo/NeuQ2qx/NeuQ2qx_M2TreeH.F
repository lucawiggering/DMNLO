************************************************

	double precision function NeuQ2qx_M2TreeH()

	implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	

	
	complex*16 DenSmZ, DenSmH(4)
	double precision colorfactor
	double precision DenSmZ1, DenSmH1(4), M2test
	double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
	integer sumU, sumUc, sumT, sumTc
	integer sumUend, sumUcend, sumTend, sumTcend

        NeuQ2qx_M2TreeH = 0d0
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


!	print *, "Tree_before=", NeuQ2qx_M2TreeH
	
! --- Mss  ---	
	
	treeSS = NC*(FR(xtype,pstype,psgen)*(dconjg(FL(xtype,pstype,psgen))*m3*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*sman +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(m1**2 - m2**2 + sman)) +
     &  Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen) +
     &  dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sman))) +
     &  dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tman) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(m1**2 + m2**2 + 2*m3**2 - tman - uman)) +
     &  AL(chimeigen,sfmeigen,sftype,sfgen)*(m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + sman*(m3**2 - uman) + m1**2*(2*m3**2 - tman - uman) + m2**2*(-2*m3**2 + sman + tman + uman))))) +
     &  FL(xtype,pstype,psgen)*(m3*dconjg(FR(xtype,pstype,psgen))*(AR(chimeigen,sfmeigen,sftype,sfgen)*(m1**2*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen)) +
     &  2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*sman + Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(-m2**2 + sman)) +
     &  AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(2*m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sman))) +
     &  dconjg(FL(xtype,pstype,psgen))*(AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(m1**2 + m3**2 - tman) +
     &  m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tman - uman)) + AR(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen)*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + sman*(m3**2 - uman) + m1**2*(2*m3**2 - tman - uman) + m2**2*(-2*m3**2 + sman + tman + uman))))))/(Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))**2

     
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeSS
!         print *, "Tree_Mss=", treeSS
!        print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH

	
        
        
! --- Mtt  ---
	do sumT = 0,sumTend
	   do sumTc = 0,sumTcend
	   
	treeTT = treeTT + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen)) +
     &	dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tman)) + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*m1*m3 +
     &  dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tman))))/((MSf(ptmeigen+sumT,pttype,ptgen)**2 - tman)*(MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tman))  
	  

           enddo
        enddo
        
     	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeTT
!    	print *, "Tree_Mtt=", treeTT
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH
     	
     	
     	
     	
! --- Muu  ---
	do sumU = 0,sumUend
	   do sumUc = 0,sumUcend
	   
	treeUU = treeUU + NC*(HR(chimeigen,pumeigen+sumU,xtype)*(dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*m1*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &	dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + 2*m3**2 + m4**2 - sman - tman) +
     &  2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*uman) + dconjg(HR(chimeigen,pumeigen+sumUc, xtype))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1 - m4)*(m1 + m4)*(m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m3**2 + m4**2 - sman))*uman +
     &  MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tman) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + uman)))) +
     &  HL(chimeigen,pumeigen+sumU, xtype)*(m1*dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) +
     &  MNeuCha(pumeigen+sumUc)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman) +
     &  2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*uman) + dconjg(HL(chimeigen,pumeigen+sumUc, xtype))*((m1 - m4)*(m1 + m4)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) + dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m3**2 + m4**2 - sman))*uman +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tman) + m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + uman)))))/((MNeuCha(pumeigen+sumU)**2 - uman)*(MNeuCha(pumeigen+sumUc)**2 - uman))



           enddo
        enddo
        
     	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeUU
!    	print *, "Tree_Muu=", treeUU
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH
     	
     	
 
     	
!--- Mst  ---

	do sumTc = 0,sumTcend

	treeST = treeST + NC*(dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(FL(xtype,pstype,psgen)*(AL(chimeigen,sfmeigen, sftype,sfgen)*Mf(pstype,psgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen)) +
     &	dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 + m3**2 - tman)) + AR(chimeigen,sfmeigen, sftype,sfgen)*(m3*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 - m2**2 + sman) + 
     &  dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*m1*(m1**2 + m2**2 + 2*m3**2 - tman - uman))) + 
     &  FR(xtype,pstype,psgen)*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*m1*m3 + 
     &  dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tman)) + AL(chimeigen,sfmeigen, sftype,sfgen)*(dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*m3*(m1**2 - m2**2 + sman) +
     &  m1*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 + m2**2 + 2*m3**2 - tman - uman)))))/((Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tman))
   
	enddo

        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeST
!     	print *, "Tree_Mst=", treeST
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH
     
       
     	
     	
! --- Mts  ---
  
	do sumT = 0,sumTend
	
	treeTS = treeTS + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*(dconjg(FL(xtype,pstype,psgen))*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(2*m1*m3*AR(chimeigen,ptmeigen+sumT,pttype,ptgen) +
     &	AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m3**2 - tman)) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m3*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 - m2**2 + sman) +
     &  AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*(m1**2 + m2**2 + 2*m3**2 - tman - uman))) +
     &  dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*m3 + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m3**2 - tman)) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m3*(m1**2 - m2**2 + sman) +
     &  m1*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m2**2 + 2*m3**2 - tman - uman)))))/((Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen+sumT,pttype,ptgen)**2 - tman))
	
	enddo
	
	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeTS
!      	print *, "Tree_Mts=", treeTS
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH
     	
   
     	
     	
! --- Msu  ---
  
	do sumUc = 0,sumUcend
	
	treeSU = treeSU + NC*(dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*(2*FL(xtype,pstype,psgen)*(AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tman) +
     &	m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + uman)) + m1*AR(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(-2*m2**2 - 2*m4**2 + sman + tman + uman))) + FR(xtype,pstype,psgen)*(2*m1*Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AL(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc)*(m1**2 - m2**2 + sman) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(2*m1**4 + (m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sman - 3*tman - uman) + (m3**2 + m4**2)*uman - uman**2 + m2**2*(-2*m4**2 + sman - tman + uman))))) +
     &  dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(2*FR(xtype,pstype,psgen)*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tman) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + uman)) + AL(chimeigen,sfmeigen,sftype,sfgen)*m1*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(-2*m2**2 - 2*m4**2 + sman + tman + uman))) + FL(xtype,pstype,psgen)*(2*AL(chimeigen,sfmeigen,sftype,sfgen)*m1*Mf(pstype,psgen)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AR(chimeigen,sfmeigen,sftype,sfgen)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sman) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(2*m1**4 + (m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sman - 3*tman - uman) +
     &  (m3**2 + m4**2)*uman - uman**2 + m2**2*(-2*m4**2 + sman - tman + uman))))))/(2.*(Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MNeuCha(pumeigen+sumUc)**2 - uman))
	
	enddo
	
	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeSU
!     	print *, "Tree_Msu=", treeSU
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH	
     	
     	
     	
     	
! --- Mus  ---
  
	do sumU = 0,sumUend
	
	treeUS = treeUS + NC*(HR(chimeigen,pumeigen+sumU,xtype)*(2*dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m3**2 - tman) + 
     &	AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(m1**2 - m4**2 + uman)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(-2*m2**2 - 2*m4**2 + sman + tman + uman))) + dconjg(FL(xtype,pstype,psgen))*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen)*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m2**2 + sman) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*m1**4 + (m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sman - 3*tman - uman) + (m3**2 + m4**2)*uman - uman**2 + m2**2*(-2*m4**2 + sman - tman + uman))))) +
     &  HL(chimeigen,pumeigen+sumU,xtype)*(2*dconjg(FL(xtype,pstype,psgen))*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m3**2 - tman) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m4**2 + uman)) + m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m2**2 + 2*m3**2 - tman - uman) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(-2*m2**2 - 2*m4**2 + sman + tman + uman))) + dconjg(FR(xtype,pstype,psgen))*(2*m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU)*(m1**2 - m2**2 + sman) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*m1**4 + (m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sman - 3*tman - uman) +
     &  (m3**2 + m4**2)*uman - uman**2 + m2**2*(-2*m4**2 + sman - tman + uman))))))/(2.*(Mf(pstype,psgen)**2 - sman + (cI*Mf(pstype,psgen)*width(pstype)))*(MNeuCha(pumeigen+sumU)**2 - uman))
	
	enddo
	
	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeUS
!      	print *, "Tree_Mus=", treeUS
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH

     	
     	
     	
! --- Mtu  ---
  
	do sumT = 0,sumTend
	   do sumUc = 0,sumUcend
	
	treeTU = treeTU + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*(dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tman) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + uman))) + dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*(m1*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tman) +
     &  m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + uman)))))/((MSf(ptmeigen+sumT,pttype,ptgen)**2 - tman)*(MNeuCha(pumeigen+sumUc)**2 - uman))
     
	   enddo
	enddo
	
	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeTU
!     	print *, "Tree_Mtu=", treeTU
!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeH	
     	     	
     	
     	
     	
! --- Mut  ---
  
	do sumTc = 0,sumTcend
	   do sumU = 0,sumUend
	
	treeUT = treeUT + NC*(dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(HR(chimeigen,pumeigen+sumU,xtype)*(dconjg(AL(chimeigen, ptmeigen+sumTc,pttype,ptgen))*m1*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen) +
     &	AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + dconjg(AR(chimeigen, ptmeigen+sumTc,pttype,ptgen))*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m3**2 - tman) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(m1**2 - m4**2 + uman))) + HL(chimeigen,pumeigen+sumU,xtype)*(m1*dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sman - tman)) + dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m3**2 - tman) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m4**2 + uman)))))/((MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tman)*(MNeuCha(pumeigen+sumU)**2 - uman))
     
	   enddo
	enddo
	
	
        NeuQ2qx_M2TreeH = NeuQ2qx_M2TreeH + treeUT
!     	print *, "Tree_Mut=", treeUT
!     	print *, "NeuQ2qx_M2TreeH=", NeuQ2qx_M2TreeH

     		
	
       end
