!**********************************************************************************************************************************

	double precision function QQ2xx_M2Tree23VH(argx,argyh)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to3.h"

	

	integer a,b,i,j,argx,argyh
	complex*16 Coupl23_A,Coupl23_At,Coupl23_Au,Coupl23_Ac,Coupl23_Atc,Coupl23_Auc,Coupl23_B,Coupl23_Bt,Coupl23_Bu,Coupl23_Bc,Coupl23_Btc,Coupl23_Buc,
     &             propa,propac,treeQbQb,treeQbS,treeQbTU,treeQbTUb,treeSS,treeTU,treeTbUb,treeTUbTU,treeSTU,treeSTUb 
        double precision sfmt,sfmu,sfmtc,sfmuc,sfwt,sfwu,sfwtc,sfwuc


        QQ2xx_M2Tree23VH = 0d0

!==================================================================================================================================	
!                          ---               Vector-Higgs-final-states                    ---
! =================================================================================================================================


!==========================================================================qb*qb===================================================
   
      do i=sumstartV, sumendV
     
        treeQbQb = 0d0
        treeQbS = 0d0
        treeQbTU = 0d0
        treeQbTUb = 0d0

        Coupl23_A = VVSfSf(ztype,i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2)
        Coupl23_B = LLR(i,argx,argyh)
        propa = DenSmV23(i)

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))
          Coupl23_Bc = dconjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV23(j))
 
        treeQbQb = treeQbQb + (4*Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc)/(propa*propac)

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + colorfactor23*dreal(treeQbQb)

!        print *, "Tree_Mqbqb=", colorfactor23*treeQbQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH

! --- MqbH+MqbV ---

!=======================================================qb*s (Higgs)+ qb*s (Vector)================================================
        

!****************************************************************Higgs*************************************************************

! - from i-combination  * - from Vectorpropagator

        do j=sumstartH, sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_Bc = dconjg(JLR(argx,j,argyh)) 
          propac = dconjg(DenSmH23(j))

!(p1(Higgs)-radiation)

	treeQbS =  treeQbS + (Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V1radc*(4*p1k1 + 8*p1k2 + sqrtS**2*(-3 - muu3**2 + muu4**2 + 3*muu5**2 + 2*x1 + x2)))/
     &  (propa*propac*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2))

!(p2(Higgs)-radiation)

	treeQbS =  treeQbS - ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V2radc*(4*p1k1 + 8*p1k2 + sqrtS**2*(3 + muu3**2 - muu4**2 - 3*muu5**2 - 4*x1 - 5*x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))))
  
!**********************************************************************************************************************************

        enddo

!*************************************************************Vector***************************************************************
!-(i-combination)  

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl23_Bc = conjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV23(j))


!(p1(Vector)-radiation)

	treeQbS =  treeQbS - ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V1radc*(8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2))) 

!(p2(Vector)-radiation)

	treeQbS =  treeQbS - ((Coupl23_A*Coupl23_Ac*Coupl23_B*Coupl23_Bc*V2radc*(8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2)))/
     &    (propa*propac*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))))
  
!**********************************************************************************************************************************

        enddo
         
        QQ2xx_M2Tree23VH= QQ2xx_M2Tree23VH + colorfactor23*dreal(treeQbS)

!    	print *, "Tree_Mqbs=", colorfactor23*treeQbS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH

! --- Mqbt+Mqbu ---

!====================================================qb*t + qb*u (sFermion)========================================================

        do j= 1,2        

            Coupl23_Atc = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc =  dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))  

!*************************************************************qb*t*****************************************************************

! - from i-combination * - from Propagator

!(p1(t)-radiation)

        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V1radc*(8*p1k1 + 12*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x2)))/
     &  (propa*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2)*
     &    (m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2 - cI*sfmtc*sfwtc))

!(p2(t)-radiation)

        treeQbTU =  treeQbTU - ((Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V2radc*(8*p1k1 + 4*p1k2 - sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2)))/
     &    (propa*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))*
     &      (-m1**2 - m3**2 + 2*p1k1 + sfmtc**2 + cI*sfmtc*sfwtc)))

!(channel(t)-radiation)

        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc*V3tradc*(8*p1k1 - 4*p1k2 + sqrtS**2*(1 - 6*muu1**2 - 2*muu2**2 - 3*muu3**2 - muu4**2 + muu5**2 + x2)))/
     &  (propa*(m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2 - cI*sfmtc*sfwtc)*
     &    (-m1**2 - m3**2 + 2*p1k1 + sfmtc**2 + cI*sfmtc*sfwtc))
  
!*************************************************************qb*u*****************************************************************
 
! - from i-combination * - from Propagator 

!(p1(u)-radiation)
     
        treeQbTU =  treeQbTU + (Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V1radc*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2)))/
     -  (propa*(2*p1k1 + 2*p1k2 + (-1 - muu1**2 + muu2**2 + muu5**2)*sqrtS**2)*
     -    (m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2 - cI*sfmuc*sfwuc))
 
!(p2(u)-radiation)
     
        treeQbTU =  treeQbTU - ((Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V2radc*(8*p1k1 + 12*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 + 4*x1 + 5*x2)))/
     -    (propa*(2*p1k1 + 2*p1k2 - sqrtS**2*(-1 + muu1**2 - muu2**2 + muu5**2 + x1 + x2))*
     -      (-m1**2 - m4**2 + 2*p1k2 + sfmuc**2 + cI*sfmuc*sfwuc)))

!(channel(u)-radiation)
     
        treeQbTU =  treeQbTU - ((Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc*V3uradc*(8*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 + x2)))/
     -    (propa*(m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2 - cI*sfmuc*sfwuc)*
     -      (-m1**2 - m4**2 + 2*p1k2 + sfmuc**2 + cI*sfmuc*sfwuc)))
  
!**********************************************************************************************************************************
   
        enddo

        QQ2xx_M2Tree23VH= QQ2xx_M2Tree23VH + colorfactor23*dreal(treeQbTU)

!    	print *, "Tree_Mqbt+Tree_Mqbu=", colorfactor23*treeQbTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH

!========================================================qb*tb + qb*ub (sFermion)==================================================
        
        do j= 1, 2  

            Coupl23_Atc = dconjg(VVSfSf(ztype,argx,j,sfmeigen1,proptypetc,sftype1,propgentc))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(VVSfSf(ztype,argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

!qb*tb
    
        treeQbTUb = treeQbTUb + (-8*Coupl23_A*Coupl23_Atc*Coupl23_B*Coupl23_Btc)/(propa*(m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2 - cI*sfmtc*sfwtc))
 
!qb*ub 

        treeQbTUb = treeQbTUb + (-8*Coupl23_A*Coupl23_Auc*Coupl23_B*Coupl23_Buc)/(propa*(m1**2 + m4**2 - 2*p1k2 - sfmuc**2 - cI*sfmuc*sfwuc))
    
!**********************************************************************************************************************************
            
        enddo

        QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + colorfactor23*dreal(treeQbTUb)

!   	print *, "Tree_MQbtb+Tree_MQbub+=",colorfactor23*treeQbTUb
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH 

      enddo

!==================================================================================================================================

! --- Mss --

!=======================================================s-channel==================================================================

!*****************************************************Higgs-exchange***************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

      treeSS = 0d0

      do i= sumstartH, sumendH   

         Coupl23_A = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
         Coupl23_B = JLR(argx,i,argyh)
         propa = DenSmH23(i)
                                
         do j= sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)) * onshell_H(j)
            Coupl23_Bc = dconjg(JLR(argx,j,argyh))
            propac = dconjg(DenSmH23(j)) 

!(p1-radiation)*(p1-radiation)

	 treeSS = treeSS + (sqrtS**2*(4*p1k1 + 4*p1k2 + (-2d0 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(-muu3**2 + 2*(-1d0 + muu4**2 + muu5**2 + x1 + x2))*V1radc*
     &                    V1rad)*Coupl23_A
     &                   *Coupl23_Ac*Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))

!(p1-radiation)*(p2-radiation)

	 treeSS = treeSS - 2*(((sqrtS**4*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-muu3**2 + 2*(-1d0 + muu4**2 + muu5**2 + x1 + x2))*V2radc*
     &                    V1rad)*Coupl23_A
     &                   *Coupl23_Ac*Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))     )

!(p2-radiation)*(p2-radiation)

	 treeSS = treeSS - ((sqrtS**2*(muu3**2 - 2*(-1d0 + muu4**2 + muu5**2 + x1 + x2))*(-4*p1k1 - 4*p1k2 + sqrtS**2*(-2d0 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    V2radc*V2rad)*Coupl23_A
     &                   *Coupl23_Ac*Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))


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
         Coupl23_B = LLR(i,argx,argyh)
         propa = DenSmV23(i)

         do j=sumstartV, sumendV

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j)
            Coupl23_Bc = dconjg(LLR(j,argx,argyh))
            propac = dconjg(DenSmV23(j)) 

!- from Vector-Vector-contraction * - from Polarisationvector 

!(p1-radiation)*(p1-radiation)
 
         treeSS = treeSS + ((4*p1k1 + 4*p1k2 + (-2d0 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 + sqrtS**2*(-1d0 + 4*muu2**2 + muu5**2 - x1 - x2))*V1radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))    

!(p1-radiation)*(p2-radiation)
 
         treeSS = treeSS - 2*((((-1d0 + 2*muu1**2 + 2*muu2**2 - muu5**2)*sqrtS**4*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))     ) 

!(p2-radiation)*(p2-radiation)
 
         treeSS = treeSS + ((4*p1k1 + 4*p1k2 - sqrtS**2*(-1d0 + 4*muu1**2 + muu5**2 + x1 + x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    V2radc*V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )) 

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
            Coupl23_Bc = dconjg(JLR(argx,j,argyh))
            propac = dconjg(DenSmH23(j)) 

!- from Vector-Vector-contraction * - from Vector-propagator * - from Polarisationvector
 
!(p1(V)-radiation)*(p1(H)-radiation)

         treeSS = treeSS - 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + 8*p1k2 + sqrtS**2*(-3 - muu3**2 + muu4**2 + 3*muu5**2 + x1 - x2))*
     &                    V1radc*V1rad)*
     &                    Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(V)-radiation)*(p2(H)-radiation)

         treeSS = treeSS + 2*((sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(-4*p1k1 - 8*p1k2 + sqrtS**2*(3 + muu3**2 - muu4**2 - 3*muu5**2 - x1 + x2))*V2radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                   (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                   (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )


!(p2(V)-radiation)*(p1(H)-radiation)

         treeSS = treeSS - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k1 - 8*p1k2 + sqrtS**2*(-muu3**2 + muu4**2 + 3*(-1 + muu5**2 + x1 + x2)))*V1radc*
     &                    V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )


!(p2(V)-radiation)*(p2(H)-radiation)

         treeSS = treeSS + 2*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*(4*p1k1 + 8*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - 3*(-1 + muu5**2 + x1 + x2)))*
     &                    V2radc*V2rad)*
     &                    Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   )


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

      QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + colorfactor23*dreal(treeSS)

!        print *, "Tree_Mss=", colorfactor23*treeSS
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH

!==================================================================================================================================

! --- Mtbtb+Mubub+Mtbub+Mubtb --

!===================================================sFermion-exchange (tb- und ub-Kanal)===========================================

      do i= 1, 2

         treeTbUb = 0d0
         treeSTUb = 0d0
         treeTUbTU = 0d0
        
         Coupl23_At = VVSfSf(ztype,argx,i,sfmeigen1,proptypet,sftype1,propgent)
         Coupl23_Bt = HSfSfu(argyh,i,sfmeigen2,proptypet,sftype2,propgent,sfgen2)
         sfmt = MSf(i,proptypet,propgent)
         sfwt = WSf(i,proptypet,propgent)

         Coupl23_Au =  HSfSfu(argyh,sfmeigen1,i,sftype1,proptypeu,sfgen1,propgenu)
         Coupl23_Bu = VVSfSf(ztype,argx,sfmeigen2,i,sftype2,proptypeu,sfgen2)
         sfmu = MSf(i,proptypeu,propgenu)
         sfwu = WSf(i,proptypeu,propgenu)

         do j=1, 2  

            Coupl23_Atc = dconjg(VVSfSf(ztype,argx,j,sfmeigen1,proptypetc,sftype1,propgentc))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(VVSfSf(ztype,argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)
             
!**********************************************************************************************************************************

!- from Polarisationvector * - from Polarisationvector 
   
!tb*tb
    
         treeTbUb = treeTbUb + (  (4*Coupl23_Btc*Coupl23_Atc*
     &                        Coupl23_Bt*Coupl23_At)/
     &                        ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )

!tb*ub

         treeTbUb = treeTbUb + 8*((Coupl23_Auc*Coupl23_Buc*
     &                        Coupl23_Bt*Coupl23_At)/
     &                        ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))))

!ub*ub 

         treeTbUb = treeTbUb + (  (4*Coupl23_Auc*Coupl23_Buc*
     &                        Coupl23_Au*Coupl23_Bu)/
     &                        ((m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*(m4**2 - 2*p1k2 + m1**2 -  (sfmu**2-cI*sfmu*sfwu)))  )

!**********************************************************************************************************************************

         enddo

         QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + colorfactor23*dreal(treeTbUb)

!   	print *, "Tree_Mtbtb+Tree_Mubub+Tree_Mubtb=", treeTbUb
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH 


!=================================================u*s (Higgs)+ u*s (Vector) + t*s (Higgs)+ t*s (Vector)============================

         do j= sumstartH, sumendH

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
            Coupl23_Bc = dconjg(JLR(argx,j,argyh))
            propac = dconjg(DenSmH23(j)) 

!**************************************************************tb*Higgs************************************************************

!- from Polarisationvector * - from Polarisationvector - from i-combination

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb - ((4*p1k1 + 8*p1k2 + sqrtS**2*(-3 - muu3**2 + muu4**2 + 3*muu5**2 + 2*x1 + x2))*
     &                          (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                          Coupl23_Bt*Coupl23_At)/
     &                          (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb - ((4*p1k1 + 8*p1k2 + sqrtS**2*(3 + muu3**2 - muu4**2 - 3*muu5**2 - 4*x1 - 5*x2))*
     &                          (-colorfactor23*((Coupl23_Ac*Coupl23_Bc*V2radc*
     &                          Coupl23_Bt*Coupl23_At)/
     &                          (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))


!**************************************************************ub*Higgs************************************************************

!- from Polarisationvector * - from Polarisationvector - from i-combination

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb - ((4*p1k1 + 8*p1k2 + sqrtS**2*(-3 - muu3**2 + muu4**2 + 3*muu5**2 + 2*x1 + x2))*
     &                          (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb - ((4*p1k1 + 8*p1k2 + sqrtS**2*(3 + muu3**2 - muu4**2 - 3*muu5**2 - 4*x1 - 5*x2))*
     &                          (-colorfactor23*((Coupl23_Ac*Coupl23_Bc*V2radc*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!**********************************************************************************************************************************

         enddo
                   
         do j= sumstartV, sumendV

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
            Coupl23_Bc = dconjg(LLR(j,argx,argyh))
            propac = dconjg(DenSmV23(j)) 

!**************************************************************tb*Vector*********************************************************** 

!- from Propagator * - from Polarisationvector * - from Polarisationvector - from i-combination  

!(p1(s)-radiation)

	 treeSTUb =  treeSTUb + (8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2))*
     &                          (colorfactor23*(Coupl23_Bc*Coupl23_Ac*V1radc*
     &                          Coupl23_Bt*Coupl23_At)/
     &                          (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

	 treeSTUb =  treeSTUb - ((8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                          (-colorfactor23*((Coupl23_Bc*Coupl23_Ac*V2radc*
     &                          Coupl23_Bt*Coupl23_At)/
     &                          (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!**************************************************************ub*Vector***********************************************************     

!- from Propagator * - from Polarisationvector * - from Polarisationvector - from i-combination                  

!(p1(s)-radiation)

          treeSTUb =  treeSTUb + (8*p1k1 + 8*p1k2 + sqrtS**2*(-4 + 2*muu1**2 + 6*muu2**2 + 2*muu5**2 - x1 - x2))*
     &                          (colorfactor23*(Coupl23_Bc*Coupl23_Ac*V1radc*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(s)-radiation)

          treeSTUb =  treeSTUb - ((8*p1k1 + 8*p1k2 - sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                          (-colorfactor23*((Coupl23_Bc*Coupl23_Ac*V2radc*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!********************************************************************************************************************************** 

         enddo

         QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + dreal(treeSTUb)

!    	print *, "Tree_Mstb+Tree_Msub=", treeSTUb
!     	print *, "Tree_M2Tree2=", QQ2xx_M2Tree23VH

!==================================================================================================================================

! --- Mtbt+Mubu+Mtbu+Mubt --

!===================================================sFermion-exchange (tb- und ub-Kanal)*(t- und u-Kanal)==========================

         do j=1, 2       

            Coupl23_Atc = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc =  dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))         
    
!******************************************************( tb*t )********************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector * - form i-combination

!(p1-radiation)
    
         treeTUbTU = treeTUbTU - ((8*p1k1 + 12*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x2))*
     &                          ( colorfactor23*(Coupl23_Btc*Coupl23_Atc*
     &                          V1radc*Coupl23_Bt*
     &                          Coupl23_At)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))

!(p2-radiation)
    
         treeTUbTU = treeTUbTU - ((-8*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &                          (-colorfactor23*((Coupl23_Btc*Coupl23_Atc*
     &                          V2radc*Coupl23_Bt*
     &                          Coupl23_At)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))) 
 
!(p3-radiation)
    
         treeTUbTU = treeTUbTU +(8*p1k1 - 4*p1k2 + sqrtS**2*(1 - 6*muu1**2 - 2*muu2**2 - 3*muu3**2 - muu4**2 + muu5**2 + x2))*
     &                          (colorfactor23*(Coupl23_Btc*Coupl23_Atc*
     &                          V3tradc*Coupl23_Bt*
     &                          Coupl23_At)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))  


 
!******************************************************( tb*u + ub*t )*************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector * - form i-combination

!(p1(u)-radiation)

         treeTUbTU = treeTUbTU - ((8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*
     &                          (colorfactor23*(Coupl23_Auc*Coupl23_Buc*
     &                          V1radc*Coupl23_Bt*
     &                          Coupl23_At)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))) 

!(p2(u)-radiation)

         treeTUbTU = treeTUbTU + (8*p1k1 + 12*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 + 4*x1 + 5*x2))*
     &                          (-colorfactor23*((Coupl23_Auc*Coupl23_Buc*
     &                          V2radc*Coupl23_Bt*
     &                          Coupl23_At)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))) 


!(p3(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((8*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 + x2))*
     &                          Coupl23_Auc*Coupl23_Buc*V3uradc*Coupl23_Bt*Coupl23_At/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!(p1(t)-radiation)

         treeTUbTU = treeTUbTU - ((8*p1k1 + 12*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x2))*
     &                          (colorfactor23*(Coupl23_Btc*Coupl23_Atc*
     &                          V1radc*Coupl23_Au*
     &                          Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))) 

!(p2(t)-radiation)

         treeTUbTU = treeTUbTU - ((-8*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &                          (-colorfactor23*((Coupl23_Btc*Coupl23_Atc*
     &                          V2radc*Coupl23_Au*
     &                          Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))) 


!(p3(t)-radiation)

         treeTUbTU = treeTUbTU + (8*p1k1 - 4*p1k2 + sqrtS**2*(1 - 6*muu1**2 - 2*muu2**2 - 3*muu3**2 - muu4**2 + muu5**2 + x2))*
     &                          (colorfactor23*(Coupl23_Btc*Coupl23_Atc*
     &                          V3tradc*Coupl23_Au*
     &                          Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                          (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!**************************************************************( ub*u )************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector * - form i-combination

!(p1(u)-radiation)

         treeTUbTU = treeTUbTU - ((8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*
     &                          (colorfactor23*(Coupl23_Auc*Coupl23_Buc*
     &                          V1radc*Coupl23_Au*
     &                          Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))) 

!(p2(u)-radiation)

         treeTUbTU = treeTUbTU + (8*p1k1 + 12*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 + 4*x1 + 5*x2))*
     &                          (-colorfactor23*((Coupl23_Auc*Coupl23_Buc*
     &                          V2radc*Coupl23_Au*
     &                          Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))) 

!(p3(u)-radiation)

         treeTUbTU = treeTUbTU - colorfactor23*((8*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 + x2))*
     &                           Coupl23_Auc*Coupl23_Buc*V3uradc*Coupl23_Au*Coupl23_Bu/
     &                           ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                           (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!**********************************************************************************************************************************

         enddo

         QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + dreal(treeTUbTU)

!   	print *, "Tree_Mtbt+Tree_Mubu+Tree_Mubt+Tree_Mtbu=", treeTUbTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH 

!==================================================================================================================================

      enddo

! --- Mtt+Muu+Mtu+Mut --

!===================================================sFermion-exchange (t- und u-Kanal)=============================================

      do i= 1, 2

         treeTU = 0d0
         treeSTU = 0d0

         Coupl23_At = VSfSft(argx,i,sfmeigen1,proptypet,sftype1,propgent,sfgen1)
         Coupl23_Bt = HSfSfu(argyh,i,sfmeigen2,proptypet,sftype2,propgent,sfgen2)
         sfmt = MSf(i,proptypet,propgent)
         sfwt = WSf(i,proptypet,propgent)

         Coupl23_Au = HSfSfu(argyh,sfmeigen1,i,sftype1,proptypeu,sfgen1,propgenu)
         Coupl23_Bu = VSfSft(argx,sfmeigen2,i,sftype2,proptypeu,sfgen2,propgenu)
         sfmu = MSf(i,proptypeu,propgenu)
         sfwu = WSf(i,proptypeu,propgenu)

         V3trad = VSfSf(ztype,i,i,proptypet,proptypet,propgent,propgent)                      !t-channel-radiation
         V3urad = VSfSf(ztype,i,i,proptypeu,proptypeu,propgenu,propgenu)                      !u-channel-radiation   

         do j= 1, 2  

            Coupl23_Atc = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc =  dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))        
  
!******************************************************( t*t )*********************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1-radiation)*(p1-radiation)
    
         treeTU = treeTU + colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + 8*p1k2 + sqrtS**2*(4*muu2**2 - muu3**2 + 2*(-1 + muu4**2 + muu5**2 - x2)))*
     &                    V1radc*V1rad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p1)-radiation)*(p2-radiation)
    
         treeTU = treeTU - 2*colorfactor23*(((sqrtS**2*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 - x2))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))     )

!(p2-radiation)*(p2-radiation)
    
         treeTU = treeTU + colorfactor23*(  ((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p1-radiation)*(p3-radiation)
    
         treeTU = treeTU + 2*colorfactor23*(  ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 - x2))*
     &                    V3tradc*V1rad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))    )

!(p2-radiation)*(p3-radiation)

    
         treeTU = treeTU - 2*colorfactor23*((((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V3tradc*
     &                    V2rad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))    )

!(p3-radiation)*(p3-radiation)
    
         treeTU = treeTU + colorfactor23*(  ((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu1**2 - 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + 2*x2))*V3tradc*
     &                    V3trad)*Coupl23_At                  
     &                   *Coupl23_Bt*Coupl23_Atc
     &                   *Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )

!**********************************************************Dipoles( t*t )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeTU - 

!p2-splitting

!         treeTU = treeTU -  

!      endif 
 
!******************************************************( t*u )*********************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1(t)-radiation)*(p1(u)-radiation)

           treeTU = treeTU + 2*colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu4**2 + muu5**2 - x1 - x2))*V1radc*V1rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &  ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)
     &    *(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))   )

!(p1(t)-radiation)*(p2(u)-radiation)

           treeTU = treeTU + 2*colorfactor23*(  -((sqrtS**4*(2*muu1**2 + 2*muu2**2 + muu3**2 - 2*muu4**2 - 2*muu5**2 - x1)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                          V2radc*V1rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &       ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &         (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &         (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &           sqrtS**2*x1 - sqrtS**2*x2)))   )

!(p2(t)-radiation)*(p1(u)-radiation)

            treeTU = treeTU + 2*colorfactor23*(   -((sqrtS**4*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                           V1radc*V2rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &       ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &         (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &         (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)
     &       *(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &           sqrtS**2*x1 - sqrtS**2*x2)))  )
 
!(p2(t)-radiation)*(p2(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(   ((4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu4**2 + muu5**2 + x1 + x2))*
     &       (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &     ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   )
 
!(p3(t)-radiation)*(p1(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(    (sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)
     &               *(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                               V1radc*V3trad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &     ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &       (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &       (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!(p1(t)-radiation)*(p3(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(          ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &       (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu4**2 + muu5**2 - x1 - x2))*V3uradc*V1rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &     ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &       (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &       (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))  )

!(p3(t)-radiation)*(p2(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(   -(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &         (4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu4**2 + muu5**2 + x1 + x2))*V2radc*V3trad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &       ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &         (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &         (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &         (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)))  )

!(p2(t)-radiation)*(p3(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(  -((sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)
     &               *(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &         V3uradc*V2rad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &       ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &         (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &         (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &         (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)))  )
 
!(p3(t)-radiation)*(p3(u)-radiation)

             treeTU = treeTU + 2*colorfactor23*(   (sqrtS**4*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(-2 + 2*muu1**2 + 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2)*
     &                           V3uradc*V3trad)*Coupl23_Buc
     &                   *Coupl23_Auc*Coupl23_At
     &                   *Coupl23_Bt/
     &     ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &       (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &       (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))  )

!**********************************************************Dipoles( t*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeTU - 

!p2-splitting

!         treeTU = treeTU -  

!      endif 
 
!**************************************************************( u*u )*************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector

!(p1(u)-radiation)*(p1(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*V1radc*
     &                    V1rad)*Coupl23_Bu               
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )


!(p1(u)-radiation)*(p2(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 + x2))*V2radc*
     &                    V1rad)*Coupl23_Bu               
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))   )

!(p2(u)-radiation)*(p2(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*(4*p1k1 + 8*p1k2 + sqrtS**2*(-4*muu1**2 + muu3**2 - 2*(-1 + muu4**2 + muu5**2 + x1 + x2)))*
     &                    V2radc*V2rad)*Coupl23_Bu              
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p1(u)-radiation)*(p3(u)-radiation)

         treeTU = treeTU + 2*colorfactor23*( ((4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*V3uradc*
     &                    V1rad)*Coupl23_Bu               
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))      )

!(p2(u)-radiation)*(p3(u)-radiation)

         treeTU = treeTU - 2*colorfactor23*( (((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 + x2))*
     &                    V3uradc*V2rad)*Coupl23_Bu              
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))    )

!(p3(u)-radiation)*(p3(u)-radiation)

         treeTU = treeTU + colorfactor23*(  ((4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*(4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu1**2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - muu5**2 - 2*x1))*
     &                    V3uradc*V3urad)*Coupl23_Bu               
     &                   *Coupl23_Au*Coupl23_Buc
     &                   *Coupl23_Auc/((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))  )

!**********************************************************Dipoles( u*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeTU - 

!p2-splitting

!         treeTU = treeTU -  

!      endif 
 
!**********************************************************************************************************************************
            
         enddo

         QQ2xx_M2Tree23VH = QQ2xx_M2Tree23VH + dreal(treeTU)

!   	print *, "Tree_Mtt+Tree_Muu+Tree_Mut=", treeTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VH 
!==================================================================================================================================

! --- Mst+Mts+Msu+Mus ---

!===============================================u*s (Higgs)+ u*s (Vector) + t*s (Higgs)+ t*s (Vector)==============================

         do j= sumstartH, sumendH

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
            Coupl23_Bc = dconjg(JLR(argx,j,argyh))
            propac = dconjg(DenSmH23(j)) 

!************************************************************t*Higgs***************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector  

!(p1(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 + 4*p1k2 + sqrtS**2*(-2 - muu3**2 + 2*muu4**2 + 2*muu5**2 + x1))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p1(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( (sqrtS**2*(2*p1k1 + 4*p1k2 + sqrtS**2*(-2 - muu3**2 + 2*muu4**2 + 2*muu5**2 + x1))*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*V2radc*
     &                        V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-2*p1k1 - 4*p1k2 + sqrtS**2*(-1 - muu4**2 + muu5**2 + x1 + x2))*V1radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p2(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*((((2*p1k1 + 4*p1k2 + sqrtS**2*(1 + muu4**2 - muu5**2 - x1 - x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2
     &                          )))   )

!(p3(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*( ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*(2*p1k1 + 4*p1k2 + sqrtS**2*(1 + muu4**2 - muu5**2 - x1 - x2))*
     &                        V1radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )
 
!(p3(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*(2*p1k1 + 4*p1k2 + sqrtS**2*(1 + muu4**2 - muu5**2 - x1 - x2))*
     &                        V2radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!**********************************************************Dipoles( t*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeSTU - 

!p2-splitting

!         treeTU = treeSTU -  

!      endif 
 
!**************************************************************u*Higgs*************************************************************

!- from Vector-Vector-contraction * - from Polarisationvector  

!(p1(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 + 4*p1k2 - sqrtS**2*(1 + muu4**2 - muu5**2 + x2))*V1radc*
     &                        V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 + 4*p1k2 - sqrtS**2*(1 + muu4**2 - muu5**2 + x2))*V2radc*
     &                        V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-2*p1k1 - 4*p1k2 + sqrtS**2*(-muu3**2 + 2*(-1 + muu4**2 + muu5**2 + x1 + x2)))*V1radc*
     &                        V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*(2*p1k1 + 4*p1k2 + sqrtS**2*(muu3**2 - 2*(-1 + muu4**2 + muu5**2 + x1 + x2)))*
     &                        V2radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )
 
!(p3(u)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*(2*p1k1 + 4*p1k2 - sqrtS**2*(1 + muu4**2 - muu5**2 + x2))*V1radc*
     &                        V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p3(u)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*(2*p1k1 + 4*p1k2 - sqrtS**2*(1 + muu4**2 - muu5**2 + x2))*V2radc*
     &                        V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!**********************************************************Dipoles( u*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeSTU - 

!p2-splitting

!         treeTU = treeSTU -  

!      endif 
 
!**********************************************************************************************************************************

         enddo
 
         do  j= sumstartV, sumendV

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
            Coupl23_Bc = dconjg(LLR(j,argx,argyh))
            propac = dconjg(DenSmV23(j)) 

!**************************************************************t*Vector************************************************************ 

!- from Propagator * - from Vector-Vector-contraction * - from Polarisationvector 

!(p1(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(8*p1k1 + 12*p1k2 + sqrtS**2*(-3 + 8*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x1 - 3*x2))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                      /(2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*( (sqrtS**2*(4*p1k2 + sqrtS**2*(-1 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - 3*muu5**2 - x1 - x2))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        V2radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))    )

!(p2(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU + 2*((sqrtS**2*(4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - muu5**2 + x1 - x2))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        V1radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU - 2*( ((8*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 8*muu1**2 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p3(t)-radiation)*(p1(s)-radiation)

	 treeSTU =  treeSTU - 2*(((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*(4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                        V1radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p3(t)-radiation)*(p2(s)-radiation)

	 treeSTU =  treeSTU + 2*(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*(8*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 8*muu1**2 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                        V2radc*V3trad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_At
     &                       *Coupl23_Bt
     &                       /(2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!**********************************************************Dipoles( t*Vectors )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeSTU - 

!p2-splitting

!         treeTU = treeSTU -  

!      endif 

!****************************************************************u*Vector**********************************************************     

!- from Propagator * - from Vector-Vector-contraction * - from Polarisationvector                 

!(p1(u)-radiation)*(p1(s)-radiation)

          treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(8*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 3*x1 - x2))*
     &                        V1radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p1(u)-radiation)*(p2(s)-radiation)

          treeSTU =  treeSTU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - muu5**2 + x1 + x2))*
     &                        V2radc*V1rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p2(u)-radiation)*(p1(s)-radiation)

          treeSTU =  treeSTU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k2 + sqrtS**2*(-1 + 4*muu1**2 + 4*muu2**2 + muu3**2 - muu4**2 - 3*muu5**2 - x1 + x2))*
     &                        V1radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(u)-radiation)*(p2(s)-radiation)

          treeSTU =  treeSTU - 2*( ((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*(8*p1k1 + 12*p1k2 - sqrtS**2*(-3 + 8*muu1**2 - muu3**2 + muu4**2 + 3*muu5**2 + 3*x1 + 3*x2))*
     &                        V2radc*V2rad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p3(u)-radiation)*(p1(s)-radiation)

          treeSTU =  treeSTU - 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*(8*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 3*x1 - x2))*
     &                        V1radc*V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p3(u)-radiation)*(p2(s)-radiation)

          treeSTU =  treeSTU + 2*( ((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*(4*p1k2 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 - muu3**2 + muu4**2 - muu5**2 + x1 + x2))*
     &                        V2radc*V3urad)*Coupl23_Ac
     &                       *Coupl23_Bc*Coupl23_Bu
     &                       *Coupl23_Au
     &                       /(2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!**********************************************************Dipoles( u*Vectors )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!         treeTU = treeSTU - 

!p2-splitting

!         treeTU = treeSTU -  

!      endif
 
!********************************************************************************************************************************** 

         enddo

         QQ2xx_M2Tree23VH= QQ2xx_M2Tree23VH + colorfactor23*dreal(treeSTU)

!    	print *, "Tree_Mst+Tree_Msu=", colorfactor23*treeSTU
!     	print *, "Tree_M2Tree2=", QQ2xx_M2Tree23VH
!==================================================================================================================================

      enddo

!       QQ2xx_M2Tree23VH = treeSTU !Nur bestimmte Kanaele sollen zaehlen
	
      end
