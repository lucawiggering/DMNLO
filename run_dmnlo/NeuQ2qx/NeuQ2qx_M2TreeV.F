************************************************

	double precision function NeuQ2qx_M2TreeV()

	implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision colorfactor
	double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
	integer sumT, sumTc, sumU, sumUc
	integer k
        

        if (xtype.eq.3) then
                k = 3
        else if (xtype.eq.4) then
                k = 1
        endif
        

        colorfactor = 1d0
        
        NeuQ2qx_M2TreeV = 0d0
        treeSS = 0d0
        treeTT = 0d0
        treeUU = 0d0
        treeST = 0d0
        treeTS = 0d0
        treeSU = 0d0
        treeUS = 0d0
        treeTU = 0d0
        treeUT = 0d0



! ==================================================================================================================================	
!                             ---                     vector boson final states           ---
! ==================================================================================================================================

	
!	print *, "Tree_before=", NeuQ2qx_M2TreeV
	
! --- Mss  ---	


        treeSS = (-1.*colorfactor*NC*(AL(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*
     &          (BR(xtype,ftype)*((m1**4*(m3**2 - 1.*sman) + m2**4*(-1.*m3**2 + sman) + sman*(m4**4 + (-1.*m3**2 + sman)*tman - 1.*m4**2*(sman + tman - 1.*uman)) + 
     &                  m1**2*(2.*m4**4 - 2.*m4**2*sman - 1.*(m3**2 - 1.*sman)*(tman + uman)) + m2**2*(-2.*m4**4 + 2.*m4**2*sman + (m3**2 - 1.*sman)*(sman + tman + uman)))*
     &                dconjg(BR(xtype,ftype)) + 3.*m3*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BL(xtype,ftype)*Mf(pstype,psgen)*(3.*m3*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(BR(xtype,ftype)) + 
     &               (m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(BL(xtype,ftype))*Mf(pstype,psgen))) + 
     &         m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(BL(xtype,ftype)*Mf(pstype,psgen)*
     &             ((2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                dconjg(BL(xtype,ftype)) + 6.*m3*m4**2*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BR(xtype,ftype)*(6.*m3*m4**2*sman*dconjg(BL(xtype,ftype)) + 
     &               (2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*dconjg(BR(xtype,ftype))*
     &                Mf(pstype,psgen)))) + AR(chimeigen,sfmeigen,sftype,sfgen)*
     &       (m1*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(BR(xtype,ftype)*Mf(pstype,psgen)*
     &             ((2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                dconjg(BR(xtype,ftype)) + 6.*m3*m4**2*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BL(xtype,ftype)*(6.*m3*m4**2*sman*dconjg(BR(xtype,ftype)) + 
     &               (2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*dconjg(BL(xtype,ftype))*
     &                Mf(pstype,psgen))) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*
     &          (BL(xtype,ftype)*((m1**4*(m3**2 - 1.*sman) + m2**4*(-1.*m3**2 + sman) + sman*(m4**4 + (-1.*m3**2 + sman)*tman - 1.*m4**2*(sman + tman - 1.*uman)) + 
     &                  m1**2*(2.*m4**4 - 2.*m4**2*sman - 1.*(m3**2 - 1.*sman)*(tman + uman)) + m2**2*(-2.*m4**4 + 2.*m4**2*sman + (m3**2 - 1.*sman)*(sman + tman + uman)))*
     &                dconjg(BL(xtype,ftype)) + 3.*m3*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BR(xtype,ftype)*Mf(pstype,psgen)*(3.*m3*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(BL(xtype,ftype)) + 
     &               (m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(BR(xtype,ftype))*Mf(pstype,psgen))))))/
     &  (m4**2*(sman - 1.*Mf(pstype,psgen)**2)**2)

        NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeSS
	
	
! --- Mtt  ---
 
        do sumT = 0, 1
          do sumTc = 0, 1
     
        treeTT = (-1.*colorfactor*NC*(m2**4 + (m4**2 - 1.*tman)**2 - 2.*m2**2*(m4**2 + tman))*CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*
     &    (-1.*AL(chimeigen,ptmeigen + sumT,pttype,fgen)*((m1**2 + m3**2 - 1.*tman)*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,fgen)) + 
     &         2.*m1*m3*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,fgen))) - 
     &      1.*AR(chimeigen,ptmeigen + sumT,pttype,fgen)*(2.*m1*m3*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,fgen)) + 
     &         (m1**2 + m3**2 - 1.*tman)*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,fgen))))*dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype)))/
     &  (m4**2*(tman - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2)*(tman - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

        
        NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeTT


          enddo
        enddo
        

! --- Muu  ---

        do sumU = 0, k
          do sumUc = 0, k
 
     
        treeUU = NC*(-1.*(AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*(dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &          (DR(chimeigen,pumeigen + sumU,sftype,xtype)*(((m1**6 - 2.*m4**6 + 2.*m4**4*sman + 2.*m4**4*tman + m1**4*(2.*m3**2 + 2.*m4**2 - 1.*sman - 1.*tman - 2.*uman) + 
     &                    2.*m4**4*uman - 2.*m4**2*sman*uman - 1.*tman*uman**2 + m3**2*(-4.*m4**4 + 2.*m4**2*uman + uman**2) - 
     &                    1.*m1**2*(m4**4 + m4**2*(sman + tman + uman) - 1.*uman*(sman + 2.*tman + uman) + m3**2*(-2.*m4**2 + 3.*uman)))*
     &                  dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 - 
     &               3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &            DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &             (-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)) + (m3*dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &               ((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              DR(chimeigen,pumeigen + sumU,sftype,xtype)*(-6.*m1*m4**2*uman*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &                 (m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2) + 
     &      AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*((m3*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &               ((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              DL(chimeigen,pumeigen + sumU,sftype,xtype)*(-6.*m1*m4**2*uman*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &                 (m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2 + 
     &         dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(DL(chimeigen,pumeigen + sumU,sftype,xtype)*
     &             (((m1**6 - 2.*m4**6 + 2.*m4**4*sman + 2.*m4**4*tman + m1**4*(2.*m3**2 + 2.*m4**2 - 1.*sman - 1.*tman - 2.*uman) + 2.*m4**4*uman - 2.*m4**2*sman*uman - 
     &                    1.*tman*uman**2 + m3**2*(-4.*m4**4 + 2.*m4**2*uman + uman**2) - 
     &                    1.*m1**2*(m4**4 + m4**2*(sman + tman + uman) - 1.*uman*(sman + 2.*tman + uman) + m3**2*(-2.*m4**2 + 3.*uman)))*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 - 
     &               3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &            DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &             (-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)))))/((-1.*uman + MNeuCha(pumeigen + sumU)**2)*(uman - 1.*MNeuCha(pumeigen + sumUc)**2))
        
     
 	NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeUU

          enddo
        enddo


     	
! --- Mst  ---

        do sumTc = 0, 1

        treeST =  (0.5*colorfactor*NC*dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype))*
     &    (AL(chimeigen,sfmeigen,sftype,sfgen)*(2.*m1*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 3.*m4**2*tman + tman**2 + m2**2*(m4**2 - 2.*tman - 1.*uman) - 1.*m4**2*uman + 
     &               tman*uman)*BR(xtype,ftype) + (-1.*m4**4 + m4**2*sman + m2**2*(-1.*m3**2 - 3.*m4**2 + sman) + m4**2*tman - 1.*sman*tman + m3**2*(-3.*m4**2 + tman) + 
     &               2.*m4**2*uman)*BL(xtype,ftype)*Mf(pstype,psgen)) - 
     &         1.*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m2**4*(-2.*m3**2 + sman + tman - 1.*uman) + 
     &               m1**2*(-1.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) + m4**2*tman + sman*tman - 1.*tman**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sman + tman - 1.*uman) - 
     &                  5.*m4**2*uman + tman*uman) + (m4**2 - 1.*tman)*(-1.*sman**2 - 1.*tman**2 + m3**2*(sman + tman - 1.*uman) + m4**2*(sman + tman - 1.*uman) + uman**2) + 
     &               m2**2*(-1.*sman**2 - 1.*sman*tman - 2.*tman**2 + 2.*m4**2*(sman - 1.*tman - 1.*uman) + m3**2*(2.*m4**2 + sman + 3.*tman - 1.*uman) + tman*uman + uman**2))*
     &             BR(xtype,ftype) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sman + m4**2*tman + m1**2*(-1.*m2**2 - 3.*m4**2 + tman) + m4**2*uman - 1.*tman*uman + m2**2*(-3.*m4**2 + uman))*
     &             BL(xtype,ftype)*Mf(pstype,psgen))) + AR(chimeigen,sfmeigen,sftype,sfgen)*
     &       (2.*m1*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*(m3*
     &             (m2**4 + 2.*m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 3.*m4**2*tman + tman**2 + m2**2*(m4**2 - 2.*tman - 1.*uman) - 1.*m4**2*uman + 
     &               tman*uman)*BL(xtype,ftype) + (-1.*m4**4 + m4**2*sman + m2**2*(-1.*m3**2 - 3.*m4**2 + sman) + m4**2*tman - 1.*sman*tman + m3**2*(-3.*m4**2 + tman) + 
     &               2.*m4**2*uman)*BR(xtype,ftype)*Mf(pstype,psgen)) - 
     &         1.*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m2**4*(-2.*m3**2 + sman + tman - 1.*uman) + 
     &               m1**2*(-1.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) + m4**2*tman + sman*tman - 1.*tman**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sman + tman - 1.*uman) - 
     &                  5.*m4**2*uman + tman*uman) + (m4**2 - 1.*tman)*(-1.*sman**2 - 1.*tman**2 + m3**2*(sman + tman - 1.*uman) + m4**2*(sman + tman - 1.*uman) + uman**2) + 
     &               m2**2*(-1.*sman**2 - 1.*sman*tman - 2.*tman**2 + 2.*m4**2*(sman - 1.*tman - 1.*uman) + m3**2*(2.*m4**2 + sman + 3.*tman - 1.*uman) + tman*uman + uman**2))*
     &             BL(xtype,ftype) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sman + m4**2*tman + m1**2*(-1.*m2**2 - 3.*m4**2 + tman) + m4**2*uman - 1.*tman*uman + m2**2*(-3.*m4**2 + uman))*
     &             BR(xtype,ftype)*Mf(pstype,psgen)))))/(m4**2*(sman - 1.*Mf(pstype,psgen)**2)*(tman - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

   	
	NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeST

        enddo



! --- Mts  ---

        do sumT = 0, 1

        treeTS = (-0.5*colorfactor*NC*CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*
     &    (dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(-2.*m1*AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 3.*m4**2*tman + tman**2 + m2**2*(m4**2 - 2.*tman - 1.*uman) - 1.*m4**2*uman + 
     &               tman*uman)*dconjg(BR(xtype,ftype)) + (-1.*m4**4 + m4**2*sman + m2**2*(-1.*m3**2 - 3.*m4**2 + sman) + m4**2*tman - 1.*sman*tman + m3**2*(-3.*m4**2 + tman) + 
     &               2.*m4**2*uman)*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &         AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m2**4*(-2.*m3**2 + sman + tman - 1.*uman) + 
     &               m1**2*(-1.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) + m4**2*tman + sman*tman - 1.*tman**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sman + tman - 1.*uman) - 
     &                  5.*m4**2*uman + tman*uman) + (m4**2 - 1.*tman)*(-1.*sman**2 - 1.*tman**2 + m3**2*(sman + tman - 1.*uman) + m4**2*(sman + tman - 1.*uman) + uman**2) + 
     &               m2**2*(-1.*sman**2 - 1.*sman*tman - 2.*tman**2 + 2.*m4**2*(sman - 1.*tman - 1.*uman) + m3**2*(2.*m4**2 + sman + 3.*tman - 1.*uman) + tman*uman + uman**2))*
     &             dconjg(BR(xtype,ftype)) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sman + m4**2*tman + m1**2*(-1.*m2**2 - 3.*m4**2 + tman) + m4**2*uman - 1.*tman*uman + 
     &               m2**2*(-3.*m4**2 + uman))*dconjg(BL(xtype,ftype))*Mf(pstype,psgen))) + 
     &      dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(-2.*m1*AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 3.*m4**2*tman + tman**2 + m2**2*(m4**2 - 2.*tman - 1.*uman) - 1.*m4**2*uman + 
     &               tman*uman)*dconjg(BL(xtype,ftype)) + (-1.*m4**4 + m4**2*sman + m2**2*(-1.*m3**2 - 3.*m4**2 + sman) + m4**2*tman - 1.*sman*tman + m3**2*(-3.*m4**2 + tman) + 
     &               2.*m4**2*uman)*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &         AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m2**4*(-2.*m3**2 + sman + tman - 1.*uman) + 
     &               m1**2*(-1.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) + m4**2*tman + sman*tman - 1.*tman**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sman + tman - 1.*uman) - 
     &                  5.*m4**2*uman + tman*uman) + (m4**2 - 1.*tman)*(-1.*sman**2 - 1.*tman**2 + m3**2*(sman + tman - 1.*uman) + m4**2*(sman + tman - 1.*uman) + uman**2) + 
     &               m2**2*(-1.*sman**2 - 1.*sman*tman - 2.*tman**2 + 2.*m4**2*(sman - 1.*tman - 1.*uman) + m3**2*(2.*m4**2 + sman + 3.*tman - 1.*uman) + tman*uman + uman**2))*
     &             dconjg(BL(xtype,ftype)) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sman + m4**2*tman + m1**2*(-1.*m2**2 - 3.*m4**2 + tman) + m4**2*uman - 1.*tman*uman + 
     &               m2**2*(-3.*m4**2 + uman))*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)))))/
     &  (m4**2*(sman - 1.*Mf(pstype,psgen)**2)*(tman - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2))

  
       NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeTS

        enddo


! --- Msu  ---

        do sumUc = 0, k

        treeSU =  NC*(0.5*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2.*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &          ((m3*BR(xtype,ftype)*(m1*(m1**4 - 4.*m4**4 + m4**2*sman + 2.*m4**2*tman + m1**2*(m2**2 + m4**2 - 1.*tman - 2.*uman) + tman*uman + uman**2 - 
     &                    1.*m2**2*(3.*m4**2 + uman))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 3.*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2 + 
     &            BL(xtype,ftype)*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)) + (dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (2.*m3*BL(xtype,ftype)*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              BR(xtype,ftype)*((-8.*m3**2*m4**2*sman - 4.*m4**4*sman + 4.*m4**2*sman**2 + 4.*m4**2*sman*tman - 1.*m1**4*(sman + tman - 1.*uman) + m3**2*sman*uman + 
     &                    m4**2*sman*uman - 1.*sman**2*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m3**2*uman**2 - 1.*m4**2*uman**2 + uman**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sman + tman) + m3**2*(8.*m4**2 - 2.*uman) + (sman + tman - 1.*uman)*uman) + 
     &                    m1**2*(-5.*m4**2*sman + sman**2 - 1.*m4**2*tman + 2.*sman*tman + tman**2 - 3.*m4**2*uman + sman*uman + tman*uman - 2.*uman**2 + 
     &                       m3**2*(-1.*sman - 3.*tman + uman) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sman - 1.*tman + uman)))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                  dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2) + 
     &      AR(chimeigen,sfmeigen,sftype,sfgen)*((dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (2.*m3*BR(xtype,ftype)*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              BL(xtype,ftype)*((-8.*m3**2*m4**2*sman - 4.*m4**4*sman + 4.*m4**2*sman**2 + 4.*m4**2*sman*tman - 1.*m1**4*(sman + tman - 1.*uman) + m3**2*sman*uman + 
     &                    m4**2*sman*uman - 1.*sman**2*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m3**2*uman**2 - 1.*m4**2*uman**2 + uman**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sman + tman) + m3**2*(8.*m4**2 - 2.*uman) + (sman + tman - 1.*uman)*uman) + 
     &                    m1**2*(-5.*m4**2*sman + sman**2 - 1.*m4**2*tman + 2.*sman*tman + tman**2 - 3.*m4**2*uman + sman*uman + tman*uman - 2.*uman**2 + 
     &                       m3**2*(-1.*sman - 3.*tman + uman) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sman - 1.*tman + uman)))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2 + 
     &         2.*dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*((m3*BL(xtype,ftype)*
     &               (m1*(m1**4 - 4.*m4**4 + m4**2*sman + 2.*m4**2*tman + m1**2*(m2**2 + m4**2 - 1.*tman - 2.*uman) + tman*uman + uman**2 - 1.*m2**2*(3.*m4**2 + uman))*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 3.*m4**2*(m1**2 - 1.*m2**2 + sman)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/
     &             m4**2 + BR(xtype,ftype)*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)))))/((sman - 1.*Mf(pstype,psgen)**2)*(uman - 1.*MNeuCha(pumeigen + sumUc)**2))
  
	NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeSU
     
        enddo



! --- Mus  ---

  	do sumU = 0, k     

        treeUS = NC*(0.5*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2.*AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &          ((m3*dconjg(BR(xtype,ftype))*(m1*(m1**4 - 4.*m4**4 + m4**2*sman + 2.*m4**2*tman + m1**2*(m2**2 + m4**2 - 1.*tman - 2.*uman) + tman*uman + uman**2 - 
     &                    1.*m2**2*(3.*m4**2 + uman))*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 3.*m4**2*(m1**2 - 1.*m2**2 + sman)*DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &            dconjg(BL(xtype,ftype))*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*DL(chimeigen,pumeigen + sumU,sftype,xtype)*
     &                  MNeuCha(pumeigen + sumU))/m4**2)) + (AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (2.*m3*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 6.*m1*m4**2*DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &              dconjg(BR(xtype,ftype))*((-8.*m3**2*m4**2*sman - 4.*m4**4*sman + 4.*m4**2*sman**2 + 4.*m4**2*sman*tman - 1.*m1**4*(sman + tman - 1.*uman) + m3**2*sman*uman + 
     &                    m4**2*sman*uman - 1.*sman**2*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m3**2*uman**2 - 1.*m4**2*uman**2 + uman**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sman + tman) + m3**2*(8.*m4**2 - 2.*uman) + (sman + tman - 1.*uman)*uman) + 
     &                    m1**2*(-5.*m4**2*sman + sman**2 - 1.*m4**2*tman + 2.*sman*tman + tman**2 - 3.*m4**2*uman + sman*uman + tman*uman - 2.*uman**2 + 
     &                       m3**2*(-1.*sman - 3.*tman + uman) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sman - 1.*tman + uman)))*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                  DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU))))/m4**2) + 
     &      dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*((AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (2.*m3*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*uman) + m4**2*uman + uman**2)*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 6.*m1*m4**2*DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &              dconjg(BL(xtype,ftype))*((-8.*m3**2*m4**2*sman - 4.*m4**4*sman + 4.*m4**2*sman**2 + 4.*m4**2*sman*tman - 1.*m1**4*(sman + tman - 1.*uman) + m3**2*sman*uman + 
     &                    m4**2*sman*uman - 1.*sman**2*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m3**2*uman**2 - 1.*m4**2*uman**2 + uman**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sman + tman) + m3**2*(8.*m4**2 - 2.*uman) + (sman + tman - 1.*uman)*uman) + 
     &                    m1**2*(-5.*m4**2*sman + sman**2 - 1.*m4**2*tman + 2.*sman*tman + tman**2 - 3.*m4**2*uman + sman*uman + tman*uman - 2.*uman**2 + 
     &                       m3**2*(-1.*sman - 3.*tman + uman) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sman - 1.*tman + uman)))*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sman) + m2**2*(m3**2 - 1.*sman) - 2.*m4**2*sman - 1.*m3**2*tman + sman*tman - 1.*m3**2*uman + sman*uman)*
     &                  DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU))))/m4**2 + 
     &         2.*AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*((m3*dconjg(BL(xtype,ftype))*
     &               (m1*(m1**4 - 4.*m4**4 + m4**2*sman + 2.*m4**2*tman + m1**2*(m2**2 + m4**2 - 1.*tman - 2.*uman) + tman*uman + uman**2 - 1.*m2**2*(3.*m4**2 + uman))*
     &                  DL(chimeigen,pumeigen + sumU,sftype,xtype) - 3.*m4**2*(m1**2 - 1.*m2**2 + sman)*DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &            dconjg(BR(xtype,ftype))*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sman - 1.*tman)*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sman) + (-1.*m3**2 + sman)*uman - 1.*m4**2*(sman - 1.*tman + uman))*DR(chimeigen,pumeigen + sumU,sftype,xtype)*
     &                  MNeuCha(pumeigen + sumU))/m4**2)))))/((sman - 1.*Mf(pstype,psgen)**2)*(uman - 1.*MNeuCha(pumeigen + sumU)**2))


!     	print *, "Tree_Mus=", treeUS
        NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeUS

        enddo

!     	print *, "Tree_M2Tree=", NeuQ2qx_M2TreeV


! --- Mtu  ---

        do sumT = 0, 1
          do sumUc = 0, k

        treeTU =  NC*(CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*(AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &       (dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(((-3.*m3**2*m4**2*sman - 1.*m4**4*sman + m4**2*sman**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tman) - 
     &                 1.*m3**2*m4**2*tman - 1.*m4**4*tman + 2.*m4**2*sman*tman + m4**2*tman**2 + m4**4*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 
     &                 1.*m4**2*uman**2 + m1**2*(3.*m4**4 - 4.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) - 5.*m4**2*tman + sman*tman + tman**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sman - 1.*tman - 1.*uman) - 2.*m4**2*uman + tman*uman) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*uman) + tman*uman + m4**2*(-1.*sman - 3.*tman + uman)))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sman))/m4**2 + sman + tman - (1.*sman*tman)/m4**2 + m3**2*(-3. + tman/m4**2) + 2.*uman)*
     &             dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &         (m3*(m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m4**2*tman + m2**2*(3.*m4**2 - 1.*uman) - 1.*m4**2*uman + tman*uman)*
     &            dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(m1*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &              dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2) + 
     &      AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m3*(m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m4**2*tman + m2**2*(3.*m4**2 - 1.*uman) - 
     &              1.*m4**2*uman + tman*uman)*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (m1*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) + dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2 + 
     &         dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(((-3.*m3**2*m4**2*sman - 1.*m4**4*sman + m4**2*sman**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tman) - 
     &                 1.*m3**2*m4**2*tman - 1.*m4**4*tman + 2.*m4**2*sman*tman + m4**2*tman**2 + m4**4*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 
     &                 1.*m4**2*uman**2 + m1**2*(3.*m4**4 - 4.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) - 5.*m4**2*tman + sman*tman + tman**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sman - 1.*tman - 1.*uman) - 2.*m4**2*uman + tman*uman) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*uman) + tman*uman + m4**2*(-1.*sman - 3.*tman + uman)))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sman))/m4**2 + sman + tman - (1.*sman*tman)/m4**2 + m3**2*(-3. + tman/m4**2) + 2.*uman)*
     &             dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))))/
     &  ((uman - 1.*MNeuCha(pumeigen + sumUc)**2)*(tman - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2))
   
	NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeTU

          enddo
        enddo



! --- Mut  ---

        do sumTc = 0, 1
          do sumU = 0, k     

        treeUT = NC*(dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype))*(dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*
     &       (AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*(((-3.*m3**2*m4**2*sman - 1.*m4**4*sman + m4**2*sman**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m3**2*m4**2*tman - 
     &                 1.*m4**4*tman + 2.*m4**2*sman*tman + m4**2*tman**2 + m4**4*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m4**2*uman**2 + 
     &                 m1**2*(3.*m4**4 - 4.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) - 5.*m4**2*tman + sman*tman + tman**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sman - 1.*tman - 1.*uman) - 2.*m4**2*uman + tman*uman) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*uman) + tman*uman + m4**2*(-1.*sman - 3.*tman + uman)))*DR(chimeigen,pumeigen + sumU,sftype,xtype))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sman))/m4**2 + sman + tman - (1.*sman*tman)/m4**2 + m3**2*(-3. + tman/m4**2) + 2.*uman)*
     &             DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &         (m3*(m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m4**2*tman + m2**2*(3.*m4**2 - 1.*uman) - 1.*m4**2*uman + tman*uman)*
     &            AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*(m1*DL(chimeigen,pumeigen + sumU,sftype,xtype) + DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2) + 
     &      dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m3*(m4**4 - 2.*m4**2*sman + m1**2*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m4**2*tman + m2**2*(3.*m4**2 - 1.*uman) - 
     &              1.*m4**2*uman + tman*uman)*AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (m1*DR(chimeigen,pumeigen + sumU,sftype,xtype) + DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &         AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*(((-3.*m3**2*m4**2*sman - 1.*m4**4*sman + m4**2*sman**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tman) - 1.*m3**2*m4**2*tman - 
     &                 1.*m4**4*tman + 2.*m4**2*sman*tman + m4**2*tman**2 + m4**4*uman + m3**2*tman*uman + m4**2*tman*uman - 1.*tman**2*uman - 1.*m4**2*uman**2 + 
     &                 m1**2*(3.*m4**4 - 4.*m4**2*sman + m3**2*(6.*m4**2 - 2.*tman) - 5.*m4**2*tman + sman*tman + tman**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sman - 1.*tman - 1.*uman) - 2.*m4**2*uman + tman*uman) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*uman) + tman*uman + m4**2*(-1.*sman - 3.*tman + uman)))*DL(chimeigen,pumeigen + sumU,sftype,xtype))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sman))/m4**2 + sman + tman - (1.*sman*tman)/m4**2 + m3**2*(-3. + tman/m4**2) + 2.*uman)*
     &             DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))))/((uman - 1.*MNeuCha(pumeigen + sumU)**2)*(tman - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

   
        NeuQ2qx_M2TreeV = NeuQ2qx_M2TreeV + treeUT

          enddo
        enddo

  
       end
