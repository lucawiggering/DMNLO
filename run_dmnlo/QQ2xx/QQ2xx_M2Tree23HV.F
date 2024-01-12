!**********************************************************************************************************************************

	double precision function QQ2xx_M2Tree23HV(argxh,argy)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to3.h"

	

	integer a,b,i,j,argxh,argy
	complex*16 Coupl23_A,Coupl23_At,Coupl23_Au,Coupl23_Ac,Coupl23_Atc,Coupl23_Auc,Coupl23_B,Coupl23_Bt,Coupl23_Bu,Coupl23_Bc,Coupl23_Btc,Coupl23_Buc,
     &             propa,propac,treeQbQb,treeQbS,treeQbTU,treeQbTUb,treeSS,treeTU,treeTUb,treeTUbTU,treeSTU,treeSTUb
        double precision sfmt,sfmu,sfmtc,sfmuc,sfwt,sfwu,sfwtc,sfwuc


        QQ2xx_M2Tree23HV = 0d0
                               
!==================================================================================================================================
!                           ---               Higgs-Vector-final-states                    ---
! =================================================================================================================================


!==========================================================================qb*qb===================================================
   
      do i=sumstartV, sumendV
     
        treeQbQb = 0d0
        treeQbS = 0d0
        treeQbTU = 0d0
        treeQbTUb = 0d0

        Coupl23_A = VVSfSf(ztype,i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2)
        Coupl23_B = LLR(i,argy,argxh) 
        propa = DenSmV23(i)

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))
          Coupl23_Bc = dconjg(LLR(j,argy,argxh))
          propac = dconjg(DenSmV23(j))
 
        treeQbQb = treeQbQb + (4*Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc)/(propa*propac)

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + colorfactor23*dreal(treeQbQb)

!        print *, "Tree_Mqbqb=", colorfactor23*treeQbQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV

! --- MqbH+MqbV ---

!=======================================================qb*s (Higgs)+ qb*s (Vector)================================================
        

!************************************************************Higgs*****************************************************************

        do j=sumstartH, sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_Bc = dconjg(JLR(argy,argxh,j)) 
          propac = dconjg(DenSmH23(j))

! - from i-combination  * - from Vectorpropagator *(additional - for terms linear in JLR)

!(p1(Higgs)-radiation)

	treeQbS =  treeQbS - (Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V1radc*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + x1 + 2*x2)))/
     &  (propa*propac*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2))

!(p2(Higgs)-radiation)

	treeQbS =  treeQbS + ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V2radc*(8*p1k1 + 4*p1k2 - sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))))
  
!**********************************************************************************************************************************

        enddo

!********************************************************Vector********************************************************************

! -(i-combination)  

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl23_Bc = dconjg(LLR(j,argy,argxh))
          propac = dconjg(DenSmV23(j))


!(p1(Vector)-radiation)

	treeQbS =  treeQbS - ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V1radc*(8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2)))

!(p2(Vector)-radiation)

	treeQbS =  treeQbS - ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V2radc*(8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))))
  
!**********************************************************************************************************************************

        enddo
         
        QQ2xx_M2Tree23HV= QQ2xx_M2Tree23HV + colorfactor23*dreal(treeQbS)

!    	print *, "Tree_Mqbs=", colorfactor23*treeQbS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV

! --- Mqbt+Mqbu ---

!====================================================qb*t + qb*u (sFermion)========================================================

        do j= 1,2        

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))

! - from i-combination * - from Propagator

!(p1(t)-radiation)

        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V1radc*(4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2)))/
     &  (propa*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2)*
     &    (m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2 - cI*sfmtc*sfwtc))

!(p2(t)-radiation)

        treeQbTU =  treeQbTU - ((Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V2radc*(12*p1k1 + 8*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2)))/
     &    (propa*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))*
     &      (-m1**2 - m3**2 + 2*p1k1 + sfmtc**2 + cI*sfmtc*sfwtc)))

!(channel(t)-radiation)

        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V3tradc*(4*p1k1 - 8*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + x1 - 4*x2)))/
     &  (propa*(m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2 - cI*sfmtc*sfwtc)*
     &    (-m1**2 - m3**2 + 2*p1k1 + sfmtc**2 + cI*sfmtc*sfwtc))
  
!***********************************************************************************************************************************
 
! - from i-combination * - from Propagator 

!(p1(u)-radiation)
     
        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V1radc*(12*p1k1 + 8*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - x1)))/
     &  (propa*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2)*
     &    (m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2 - cI*sfmuc*sfwuc))
 
!(p2(u)-radiation)
     
        treeQbTU =  treeQbTU - ((Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V2radc*(4*p1k1 + 8*p1k2 - sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2)))/
     &    (propa*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))*
     &      (-m1**2 - m4**2 + 2*p1k2 + sfmuc**2 + cI*sfmuc*sfwuc)))

!(channel(u)-radiation)
     
        treeQbTU =  treeQbTU  - ((Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V3uradc*(4*p1k1 - 8*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - x1)))/
     &    (propa*(m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2 - cI*sfmuc*sfwuc)*
     &      (-m1**2 - m4**2 + 2*p1k2 + sfmuc**2 + cI*sfmuc*sfwuc))) 
  
!***********************************************************************************************************************************
   
        enddo

        QQ2xx_M2Tree23HV= QQ2xx_M2Tree23HV + colorfactor23*dreal(treeQbTU)

!    	print *, "Tree_Mqbt+Tree_Mqbu=", colorfactor23*treeQbTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV

!========================================================qb*tb + qb*ub (sFermion)====================================================
        
        do j= 1, 2  

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(VVSfSf(ztype,argy,sfmeigen2,j,sftype2,proptypetc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(VVSfSf(ztype,argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)
   
!qb*tb
    
        treeQbTUb = treeQbTUb + (-8*Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc)/(propa*(m1**2 + m3**2 - 2*p1k1 - sfmtc**2 - cI*sfmtc*sfwtc))
 
!qb*ub 

        treeQbTUb = treeQbTUb + (-8*Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc)/(propa*(m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2 - cI*sfmuc*sfwuc))
    
!**********************************************************************************************************************************
            
        enddo

        QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + colorfactor23*dreal(treeQbTUb)

!   	print *, "Tree_MQbtb+Tree_MQbub+=",colorfactor23*treeQbTUb
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV 

      enddo

!=======================================================s-channel==================================================================

!*****************************************************Higgs-exchange***************************************************************

      treeSS = 0d0

      do i= sumstartH, sumendH   

         Coupl23_A = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
         Coupl23_B = JLR(argy,argxh,i)
         propa = DenSmH23(i)
                                 
         do j= sumstartH, sumendH  

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))  * onshell_H(j)
            Coupl23_Bc = dconjg(JLR(argy,argxh,j))
            propac = dconjg(DenSmH23(j)) 

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1-radiation)*(p1-radiation)

	treeSS = treeSS + (sqrtS**2*(4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*muu3**2 - muu4**2 + 2*(-1 + muu5**2 + x1 + x2))*
     &                    V1radc*V1rad)*Coupl23_A
     &                   *Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!(p1-radiation)*(p2-radiation)

	treeSS = treeSS - 2*(((sqrtS**4*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(2*muu3**2 - muu4**2 + 2*(-1 + muu5**2 + x1 + x2))*
     &                    V2radc*V1rad)*Coupl23_A
     &                   *Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))     )
 
!(p2-radiation)*(p2-radiation)

	treeSS = treeSS + (sqrtS**2*(2*muu3**2 - muu4**2 + 2*(-1 + muu5**2 + x1 + x2))*(-4*p1k1 - 4*p1k2 + sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    V2radc*V2rad)*Coupl23_A
     &                   *Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))


!**********************************************************Dipoles( Higgs*Higgs )**************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treeSS - 

!p2-splitting

!        treeSS = treeSS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
      enddo

!****************************************************Vector-exchange***************************************************************
    
      do i=sumstartV, sumendV 

         Coupl23_A = VSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
         Coupl23_B = LLR(i,argy,argxh)
         propa = DenSmV23(i)

         do j=sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j)
            Coupl23_Bc = dconjg(LLR(j,argy,argxh))
            propac = dconjg(DenSmV23(j)) 

!- from Vector-Vector-contraction * - from Polarisationvector  

!(p1-radiation)*(p1-radiation)
 
        treeSS = treeSS +  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                     (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 + muu5**2 - x1 - x2))*V1radc*
     &                     V1rad)*Coupl23_A*Coupl23_Ac
     &                    *Coupl23_B*Coupl23_Bc/
     &                     (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!(p1-radiation)*(p2-radiation)

        treeSS = treeSS - 2*((((-1d0 + 2*muu1**2 + 2*muu2**2 - muu5**2)*sqrtS**4*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                    V2radc*V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))     )


!(p2-radiation)*(p2-radiation)
 
        treeSS = treeSS + ((4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 + muu5**2 + x1 + x2))*
     &                    (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                    sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))

!**********************************************************Dipoles( Vector*Vector )************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treeSS - 

!p2-splitting

!        treeSS = treeSS -  

!      endif
  
!**********************************************************************************************************************************
        
        enddo   

!********************************************************Vector-Higgs *************************************************************

        do j=sumstartH, sumendH 

           Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
           Coupl23_Bc = dconjg(JLR(argy,argxh,j))
           propac = dconjg(DenSmH23(j)) 

!- from Vector-Vector-contraction * - from Vector-propagator * - from Polarisationvector *(additional - for terms linear in JLR)
 
!(p1(V)-radiation)*(p1(H)-radiation)

        treeSS = treeSS + 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                    (8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 - x1 + x2))*V1radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p1(V)-radiation)*(p2(H)-radiation)

        treeSS = treeSS - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 - x1 + x2))*
     &                    V2radc*V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                   (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))       )


!(p2(V)-radiation)*(p1(H)-radiation)

        treeSS = treeSS + 2*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-8*p1k1 - 4*p1k2 + sqrtS**2*(muu3**2 - muu4**2 + 3*(-1 + muu5**2 + x1 + x2)))*
     &                    V1radc*V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))      )


!(p2(V)-radiation)*(p2(H)-radiation)

        treeSS = treeSS - 2*( ((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    (8*p1k1 + 4*p1k2 - sqrtS**2*(muu3**2 - muu4**2 + 3*(-1 + muu5**2 + x1 + x2)))*V2radc*
     &                    V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                    sqrtS**2*x2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))      )


!**********************************************************Dipoles( Vector*Higgs )*************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treeSS - 

!p2-splitting

!        treeSS = treeSS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
      enddo

!**********************************************************************************************************************************

      QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + colorfactor23*dreal(treeSS)

!        print *, "Tree_Mss=", treeSS
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV

!===================================================sFermion-exchange (tb- und ub-Kanal)===========================================

      do i= 1, 2

         treeTUb = 0d0
         treeTUbTU = 0d0
         treeSTUb = 0d0
        
         Coupl23_At = HSfSft(argxh,sfmeigen1,i,sftype1,proptypet,sfgen1,propgent)
         Coupl23_Bt = VVSfSf(ztype,argy,sfmeigen2,i,sftype2,proptypet,sfgen2)
         sfmt = MSf(i,proptypet,propgent)
         sfwt = WSf(i,proptypet,propgent)

         Coupl23_Au = VVSfSf(ztype,argy,i,sfmeigen1,proptypeu,sftype1,propgenu)
         Coupl23_Bu = HSfSft(argxh,i,sfmeigen2,proptypeu,sftype2,propgenu,sfgen2)
         sfmu = MSf(i,proptypeu,propgenu)
         sfwu = WSf(i,proptypeu,propgenu)

         do j= 1, 2  

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(VVSfSf(ztype,argy,sfmeigen2,j,sftype2,proptypetc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(VVSfSf(ztype,argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

!- from Vector-Vector-contraction * - from Polarisationvector
   
!tb*tb
    
        treeTUb = treeTUb + (  (4*Coupl23_Atc*
     &                      Coupl23_Btc*Coupl23_At*
     &                      Coupl23_Bt)/
     &                      ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )
 
!tb*ub + ub*tb

        treeTUb = treeTUb + 8*((Coupl23_Buc*
     &                      Coupl23_Auc*Coupl23_At*
     &                      Coupl23_Bt)/
     &                      ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))
  
!ub*ub

        treeTUb = treeTUb + (  (4*Coupl23_Buc*
     &                      Coupl23_Auc*Coupl23_Bu*
     &                      Coupl23_Au)/
     &                      ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))  )
  
!**********************************************************************************************************************************
            
        enddo

        QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + colorfactor23*dreal(treeTUb)

!   	print *, "Tree_Mtbtb+Tree_Mubub+Tree_Mubtb=", treeTUb
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV 

!==================================================================================================================================


! --- Mstb+Mubs+Msub+Mtbs ---

!==========================================ub*s (Higgs)+ ub*s (Vector) + tb*s (Higgs)+ tb*s (Vector)===============================

         do j= sumstartH, sumendH 

!************************************************************tb*Higgs**************************************************************

           Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
           Coupl23_Bc = dconjg(JLR(argy,argxh,j))
           propac = dconjg(DenSmH23(j)) 

!(- from Vector-Vector-contraction)^2 * - from i-factors * - for terms linear in JLR   

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + x1 + 2*x2))*
     &                          ((Coupl23_Ac*Coupl23_Bc*
     &                          V1radc*Coupl23_At*
     &                          Coupl23_Bt)/
     &                          (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 4*p1k2 - sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &                          (-((Coupl23_Ac*Coupl23_Bc*
     &                          V2radc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))) )

!**********************************************************Dipoles( tb*Higgs )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTUb = treeSTUb - 

!p2-splitting

!        treeSTUb = treeSTUb -  

!      endif
  
!**************************************************************ub*Higgs************************************************************

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + x1 + 2*x2))*
     &                          ((Coupl23_Ac*Coupl23_Bc*
     &                          V1radc*Coupl23_Bu*
     &                          Coupl23_Au)/
     &                          (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 4*p1k2 - sqrtS**2*(-3 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &                          (-((Coupl23_Ac*Coupl23_Bc*
     &                          V2radc*
     &                          Coupl23_Bu*Coupl23_Au)/
     &                          (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))) )

!**********************************************************Dipoles( ub*Higgs )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTUb = treeSTUb - 

!p2-splitting

!        treeSTUb = treeSTUb -  

!      endif
  
!**********************************************************************************************************************************

         enddo 

!**************************************************************tb*Vector***********************************************************
   
         do  j= sumstartV, sumendV 

             Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
             Coupl23_Bc = dconjg(LLR(j,argy,argxh))
             propac = dconjg(DenSmV23(j)) 

!- from Propagator * (- from Vector-Vector-contraction)^2 * - from i-factors 

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2))*
     &                          ((Coupl23_Bc*Coupl23_Ac*
     &                          V1radc*Coupl23_At*
     &                          Coupl23_Bt)/
     &                          (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb - colorfactor23*((8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                          (-((Coupl23_Bc*Coupl23_Ac*
     &                          V2radc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)   ))))

!**********************************************************Dipoles( tb*Vectors )***************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTUb = treeSTUb - 

!p2-splitting

!        treeSTUb = treeSTUb -  

!      endif
  
!****************************************************************ub*Vector*********************************************************                     

!(p1(s)-radiation)

         treeSTUb =  treeSTUb + colorfactor23*(8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2))*
     &                          ((Coupl23_Bc*Coupl23_Ac*
     &                          V1radc*Coupl23_Bu*
     &                          Coupl23_Au)/
     &                          (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

         treeSTUb =  treeSTUb - colorfactor23*((8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                          (-((Coupl23_Bc*Coupl23_Ac*
     &                          V2radc*
     &                          Coupl23_Bu*Coupl23_Au)/
     &                          (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)  ))))

!**********************************************************Dipoles( ub*Vectors )***************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTUb = treeSTUb - 

!p2-splitting

!        treeSTUb = treeSTUb -  

!      endif
 
!**********************************************************************************************************************************

        enddo

        QQ2xx_M2Tree23HV= QQ2xx_M2Tree23HV + dreal(treeSTUb)

!    	print *, "Tree_Mstb+Tree_Msub=", treeSTUb
!     	print *, "Tree_M2Tree23HV=", QQ2xx_M2Tree23HV

!==================================================================================================================================


! --- Mtbt+Mubu+Mtbu+Mubt --

!===================================================sFermion-exchange (tb- und ub-Kanal)*(t- und u-Kanal)==========================

!******************************************************( tb*t )********************************************************************

         do j=1, 2   

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))

!(- from Vector-Vector-contraction)^2 * - from i-factors
   
!(p1-radiation)
     
         treeTUbTU = treeTUbTU - colorfactor23*((4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*
     &                          ((Coupl23_Atc*
     &                          V1radc*
     &                          Coupl23_Btc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2-radiation)
    
         treeTUbTU = treeTUbTU + colorfactor23*(12*p1k1 + 8*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &                          (-((Coupl23_Atc*
     &                          V2radc*
     &                          Coupl23_Btc*Coupl23_At*
     &                          Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))  )
 
!(p3-radiation)
    
         treeTUbTU = treeTUbTU + colorfactor23*(4*p1k1 - 8*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + x1 - 4*x2))*
     &                          ((Coupl23_Atc*
     &                          V3tradc*
     &                          Coupl23_Btc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!**********************************************************Dipoles( tb*t )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUbTU = treeTUbTU - 

!p2-splitting

!        treeTUbTU = treeTUbTU -  

!      endif
  
!******************************************************( tb*u  )*******************************************************************

!(p1(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((12*p1k1 + 8*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - x1))*
     &                          ((Coupl23_Buc*
     &                          V1radc*
     &                          Coupl23_Auc*Coupl23_At*
     &                          Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((-4*p1k1 - 8*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &                          (-((Coupl23_Buc*
     &                          V2radc*
     &                          Coupl23_Auc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)  ))))

!(p3(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((4*p1k1 - 8*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - x1))*
     &                          ((Coupl23_Buc*
     &                          V3uradc*
     &                          Coupl23_Auc*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))))

!**********************************************************Dipoles( tb*u )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUbTU = treeTUbTU - 

!p2-splitting

!        treeTUbTU = treeTUbTU -  

!      endif

!******************************************************( ub*t  )*******************************************************************

!(p1(t)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*
     &                          ((Coupl23_Atc*
     &                          V1radc*
     &                          Coupl23_Btc*
     &                          Coupl23_Bu*Coupl23_Au)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2(t)-radiation)

         treeTUbTU = treeTUbTU + colorfactor23* (12*p1k1 + 8*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &                          (-((Coupl23_Atc*
     &                          V2radc*
     &                          Coupl23_Btc*Coupl23_Bu*
     &                          Coupl23_Au)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))) )

!(p3(t)-radiation)

         treeTUbTU = treeTUbTU + colorfactor23*(4*p1k1 - 8*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + x1 - 4*x2))*
     &                          ((Coupl23_Atc*
     &                          V3tradc*
     &                          Coupl23_Btc*
     &                          Coupl23_Bu*Coupl23_Au)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!**********************************************************Dipoles( ub*t )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUbTU = treeTUbTU - 

!p2-splitting

!        treeTUbTU = treeTUbTU -  

!      endif
    
**************************************************************( ub*u )************************************************************

!(p1(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((12*p1k1 + 8*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - x1))*
     &                          ((Coupl23_Buc*
     &                          V1radc*
     &                          Coupl23_Auc*Coupl23_Bu*
     &                          Coupl23_Au)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))


!(p2(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((-4*p1k1 - 8*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &                          (-((Coupl23_Buc*
     &                          V2radc*
     &                          Coupl23_Auc*
     &                          Coupl23_Bu*Coupl23_Au)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)  ))))

!(p3(u)-radiation)

          treeTUbTU = treeTUbTU - colorfactor23*((4*p1k1 - 8*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - x1))*
     &                          Coupl23_Auc*
     &                          V3tradc*
     &                          Coupl23_Buc*
     &                          Coupl23_Au*Coupl23_Bu/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                         (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))


!**********************************************************Dipoles( ub*u )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUbTU = treeTUbTU - 

!p2-splitting

!        treeTUbTU = treeTUbTU -  

!      endif
  
!**********************************************************************************************************************************
            
         enddo

         QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + dreal(treeTUbTU)

!   	print *, "Tree_Mtbt+Tree_Mubu+Tree_Mubt+Tree_Mtbu=", treeTUbTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV 

!==================================================================================================================================

      enddo

! --- Mtt+Muu+Mtu+Mut --

!===================================================sFermion-exchange (t- und u-Kanal)=============================================

      do i= 1, 2

         treeTU = 0d0
         treeSTU = 0d0

         Coupl23_At = HSfSft(argxh,sfmeigen1,i,sftype1,proptypet,sfgen1,propgent)
         Coupl23_Bt = VSfSfu(argy,sfmeigen2,i,sftype2,proptypet,sfgen2,propgent)
         sfmt = MSf(i,proptypet,propgent)
         sfwt = WSf(i,proptypet,propgent)

         Coupl23_Au = VSfSfu(argy,i,sfmeigen1,proptypeu,sftype1,propgenu,sfgen1)
         Coupl23_Bu = HSfSft(argxh,i,sfmeigen2,proptypeu,sftype2,propgenu,sfgen2)
         sfmu = MSf(i,proptypeu,propgenu)
         sfwu = WSf(i,proptypeu,propgenu)

         V3trad = VSfSf(ztype,i,i,proptypet,proptypet,propgent,propgent)                      !t-channel-radiation
         V3urad = VSfSf(ztype,i,i,proptypeu,proptypeu,propgenu,propgenu)                      !u-channel-radiation   

         do j= 1, 2  

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))
   
!******************************************************( t*t )*********************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1-radiation)*(p1-radiation)
    
         treeTU = treeTU + colorfactor23*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*
     &                    V1radc*V1rad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )) )

!(p1-radiation)*(p2-radiation)
    
         treeTU = treeTU + 2*colorfactor23*( (sqrtS**2*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 + x1))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                    V2radc*V1rad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  ))    )

!(p2-radiation)*(p2-radiation)
    
         treeTU = treeTU + colorfactor23*(  ((8*p1k1 + 4*p1k2 + sqrtS**2*(2 - 4*muu1**2 - 2*muu3**2 + muu4**2 - 2*muu5**2 - 2*x1 - 2*x2))*
     &                    (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(p1-radiation)*(channel-radiation)
    
         treeTU = treeTU + 2*colorfactor23*( ((4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                    V3tradc*V1rad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p2-radiation)*(channel-radiation)

    
         treeTU = treeTU - 2*colorfactor23*((((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                    (4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 + x1))*V3tradc*
     &                    V2rad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))       )

!(channel-radiation)*(channel-radiation)
    
         treeTU = treeTU + colorfactor23*(  ((4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*
     &                    (-4*p1k1 + 4*p1k2 + sqrtS**2*(2*muu1**2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - muu5**2 - 2*x2))*
     &                    V3tradc*V3trad)*Coupl23_Bt                  
     &                   *Coupl23_At*Coupl23_Btc
     &                   *Coupl23_Atc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )

!**********************************************************Dipoles( t*t )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU -  

!      endif
  
!******************************************************( t*u )*********************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1(t)-radiation)*(p1(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu3**2 + muu5**2 - x1 - x2))*V1radc*
     &                    V1rad)* Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))             )

!(p1(t)-radiation)*(p2(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( ((sqrtS**4*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &                    V2radc*V1rad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))              )

!(p2(t)-radiation)*(p1(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( ((sqrtS**4*(2*muu1**2 + 2*muu2**2 - 2*muu3**2 + muu4**2 - 2*muu5**2 - x2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                    V1radc*V2rad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))             )

!(p2(t)-radiation)*(p2(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( ((4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu3**2 + muu5**2 + x1 + x2))*
     &                    (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))              )

!(channel(t)-radiation)*(p1(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu3**2 + muu5**2 - x1 - x2))*V1radc*
     &                    V3trad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))              )

!(p1(t)-radiation)*(channel(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( (sqrtS**2*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &                    V3uradc*V1rad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))              )

!(channel(t)-radiation)*(p2(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*(  (sqrtS**2*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &                    V2radc*V3trad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))             )

!(p2(t)-radiation)*(channel(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( (((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                    (4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu3**2 + muu5**2 + x1 + x2))*V3uradc*
     &                    V2rad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))              )


!(channel(t)-radiation)*(channel(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( (sqrtS**4*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*(-2 + 2*muu1**2 + 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2)*
     &                    V3uradc*V3trad)*Coupl23_Bt
     &                   *Coupl23_At*Coupl23_Auc
     &                   *Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))              )

!**********************************************************Dipoles( t*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU -  

!      endif
  
!**************************************************************( u*u )*************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1(u)-radiation)*(p1(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                    (8*p1k1 + 4*p1k2 + sqrtS**2*(-2 + 4*muu2**2 + 2*muu3**2 - muu4**2 + 2*muu5**2 - 2*x1))*
     &                    V1radc*V1rad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )


!(p1(u)-radiation)*(p2(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( ((sqrtS**2*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 - x1))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                    V2radc*V1rad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))   )

!(p2(u)-radiation)*(p2(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    V2radc*V2rad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(p1(u)-radiation)*(channel(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                    (4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 - x1))*V3uradc*
     &                    V1rad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(p2(u)-radiation)*(channel(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( (((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                    V3uradc*V2rad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                    (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))    )

!(channel(u)-radiation)*(channel(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(-4*p1k1 + 4*p1k2 + 
     &                    sqrtS**2*(-2*muu1**2 - 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + 2*x1))*V3uradc*
     &                    V3urad)*Coupl23_Au                        
     &                   *Coupl23_Bu*Coupl23_Auc
     &                   *Coupl23_Buc/
     &                    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))  )

!**********************************************************Dipoles( u*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU -  

!      endif
  
!**********************************************************************************************************************************
            
         enddo

         QQ2xx_M2Tree23HV = QQ2xx_M2Tree23HV + dreal(treeTU)

!   	print *, "Tree_Mtt+Tree_Muu+Tree_Mut=", treeTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HV 
!==================================================================================================================================

! --- Mst+Mts+Msu+Mus ---

!==============================================u*s (Higgs)+ u*s (Vector) + t*s (Higgs)+ t*s (Vector)===============================

!************************************************************t*Higgs***************************************************************

         do j= sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
            Coupl23_Bc = dconjg(JLR(argy,argxh,j))
            propac = dconjg(DenSmH23(j)) 

!- from Vector-Vector-contraction * - from Polarisationvector (*-1 for terms linear in JLR) 

!(p1(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + 2*p1k2 - sqrtS**2*(1 + muu3**2 - muu5**2 + x1))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p1(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*( (sqrtS**2*(4*p1k1 + 2*p1k2 - sqrtS**2*(1 + muu3**2 - muu5**2 + x1))*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*
     &                        V2radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))            )

!(p2(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        (-4*p1k1 - 2*p1k2 + sqrtS**2*(2*muu3**2 - muu4**2 + 2*(-1 + muu5**2 + x1 + x2)))*V1radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))           )

!(p2(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        (4*p1k1 + 2*p1k2 + sqrtS**2*(-2*muu3**2 + muu4**2 - 2*(-1 + muu5**2 + x1 + x2)))*V2radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))          )

!(channel(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + 2*p1k2 - sqrtS**2*(1 + muu3**2 - muu5**2 + x1))*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                        V1radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))            )
 
!(channel(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*(  ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*(4*p1k1 + 2*p1k2 - sqrtS**2*(1 + muu3**2 - muu5**2 + x1))*
     &                        V2radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))           )

!**********************************************************Dipoles( t*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTU = treeSTU - 

!p2-splitting

!        treeSTU = treeSTU -  

!      endif
  
!**************************************************************u*Higgs*************************************************************

!(p1(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                        (4*p1k1 + 2*p1k2 + sqrtS**2*(-2 + 2*muu3**2 - muu4**2 + 2*muu5**2 + x2))*V1radc*
     &                        V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))        )

!(p1(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(4*p1k1 + 2*p1k2 + sqrtS**2*(-2 + 2*muu3**2 - muu4**2 + 2*muu5**2 + x2))*
     &                        V2radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  )))        )

!(p2(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k1 - 2*p1k2 + sqrtS**2*(-1 - muu3**2 + muu5**2 + x1 + x2))*
     &                        V1radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))        )

!(p2(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*((((4*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu3**2 - muu5**2 - x1 - x2))*
     &                        (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu/                   
     &                        (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))         )

!(channel(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                        (4*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu3**2 - muu5**2 - x1 - x2))*V1radc*
     &                        V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))        )

!(channel(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*(((4*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu3**2 - muu5**2 - x1 - x2))*
     &                        (4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V2radc*
     &                        V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))         )

!**********************************************************Dipoles( u*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTU = treeSTU - 

!p2-splitting

!        treeSTU = treeSTU -  

!      endif
  
!**********************************************************************************************************************************

         enddo
   
         do j= sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
            Coupl23_Bc = dconjg(LLR(j,argy,argxh))
            propac = dconjg(DenSmV23(j)) 

!**************************************************************t*Vector************************************************************

!- from Propagator * - from Vector-Vector-contraction * - from Polarisationvector 

!(p1(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                        (4*p1k1 + 8*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 3*x2))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(p1(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        (-4*p1k1 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - muu5**2 + x1 + x2))*
     &                        V2radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))     )

!(p2(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*((sqrtS**2*(-4*p1k1 + sqrtS**2*(-1 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - 3*muu5**2 + x1 - x2))*
     &                        (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*V2rad
     &                        )*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))     )

!(p2(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        (12*p1k1 + 8*p1k2 - sqrtS**2*(-3 + 8*muu1**2 + muu3**2 - muu4**2 + 3*muu5**2 + 3*x1 + 3*x2))*
     &                        V2radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!(channel(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*(((4*p1k1 + 8*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 3*x2))*
     &                        (4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*V1radc*
     &                        V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(channel(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                        (4*p1k1 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - muu5**2 + x1 + x2))*
     &                        V2radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bt
     &                       *Coupl23_At
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )
 
!**********************************************************Dipoles( t*Vectors )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTU = treeSTU - 

!p2-splitting

!        treeSTU = treeSTU -  

!      endif
  
!****************************************************************u*Vector**********************************************************                   

!(p1(u)-radiation)*(p1(s)-radiation)

         treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                        (12*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 8*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - 3*x1 - x2))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p1(u)-radiation)*(p2(s)-radiation)

         treeSTU =  treeSTU + 2*((sqrtS**2*(4*p1k1 + sqrtS**2*(-1 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - 3*muu5**2 - x1 - x2))*
     &                        (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad
     &                        )*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))       )

!(p2(u)-radiation)*(p1(s)-radiation)

         treeSTU =  treeSTU + 2*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        (4*p1k1 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - muu5**2 - x1 + x2))*
     &                        V1radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))      )

!(p2(u)-radiation)*(p2(s)-radiation)

         treeSTU =  treeSTU - 2*(((4*p1k1 + 8*p1k2 - sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*
     &                        (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))       )

!(channel(u)-radiation)*(p1(s)-radiation)

         treeSTU =  treeSTU - 2*( ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                        (4*p1k1 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - muu5**2 - x1 + x2))*
     &                        V1radc*V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))       )

!(channel(u)-radiation)*(p2(s)-radiation)

         treeSTU =  treeSTU + 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                        (4*p1k1 + 8*p1k2 - sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*
     &                        V2radc*V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Au
     &                       *Coupl23_Bu
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))       )

!**********************************************************Dipoles( u*Vectors )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSTU = treeSTU - 

!p2-splitting

!        treeSTU = treeSTU -  

!      endif
 
!**********************************************************************************************************************************

         enddo

         QQ2xx_M2Tree23HV= QQ2xx_M2Tree23HV + colorfactor23*dreal(treeSTU)

!    	print *, "Tree_Mst+Tree_Msu=", treeSTU
!     	print *, "Tree_M2Tree23HV=", QQ2xx_M2Tree23HV

!==================================================================================================================================

      enddo   

!       QQ2xx_M2Tree23HV = treeSTU   !Nur bestimmte Kanaele sollen zaehlen

   	
      end
