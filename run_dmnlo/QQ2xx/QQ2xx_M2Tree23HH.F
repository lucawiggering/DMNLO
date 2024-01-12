!**********************************************************************************************************************************

	double precision function QQ2xx_M2Tree23HH(argxh,argyh)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to3.h"

	

	integer a,b,c,d,i,j,argxh,argyh
	complex*16 treeQQ,treeQbQb,treeQbQ,treeSS,treeQS,treeQbS,treeTUQ,treeQbTU,treeTUS,treeTU,Coupl23_A,
     &             Coupl23_At,Coupl23_Au,Coupl23_Ac,Coupl23_Atc,Coupl23_Auc,Coupl23_B,Coupl23_Bt,Coupl23_Bu,Coupl23_Bc,Coupl23_Btc,Coupl23_Buc,
     &             propa,propac
	double precision sfmt,sfmu,sfmtc,sfmuc,sfwt,sfwu,sfwtc,sfwuc


        QQ2xx_M2Tree23HH = 0d0                        
 
! =================================================================================================================================	
!                             ---                     Higgs-Higgs-final-states                    ---
! =================================================================================================================================



!===============================================quartic sfermion-sfermion-Higgs-Higgs==============================================
!- from Polarisationvector

        treeQQ = (0d0,0d0)
        treeQS = (0d0,0d0)

        Coupl23_A = HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)
        Coupl23_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))
    
!(p1-radiation)*(p1-radiation)

        treeQQ = treeQQ - ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*V1rad)
     &                   *Coupl23_A*Coupl23_Ac/
     &                    ((2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))

!(p1-radiation)*(p2-radiation)+(p2-radiation)*(p1-radiation)

        treeQQ = treeQQ - 2d0*( -((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad)
     &                   *Coupl23_A*Coupl23_Ac/
     &                    ((2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p2-radiation)*(p2-radiation)

        treeQQ = treeQQ + (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)    
     &                   *Coupl23_A*Coupl23_Ac/
     &                    ((2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))

!**********************************************************Dipoles( q*q )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeQQ = treQQ - 

!p2-splitting

!        treeQQ = treQQ -  

!     endif
  
!**********************************************************************************************************************************

      QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQQ)

!        print *, "Tree_Mqq=", treeQQ
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!=======================================================q*s (Higgs)+ q*s (Vector)==================================================

!- from Polarisationvector * - from i-combination
        
!Higgs

         do j= sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))  
            Coupl23_Bc = dconjg(ILR(j,argxh,argyh))
            propac = dconjg(DenSmH23(j)) 

!(p1(Q)-radiation)*(p1(H)-radiation)

	 treeQS =  treeQS + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*V1rad)
     &                     *Coupl23_A* Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(Q)-radiation)*(p2(H)-radiation)+(p2(H)-radiation)*(p1(Q)-radiation)

	 treeQS =  treeQS - 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad)
     &                     *Coupl23_A* Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))) )

!(p2(Q)-radiation)*(p1(H)-radiation)+(p1(H)-radiation)*(p2(Q)-radiation)

	 treeQS =  treeQS - 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*V2rad)
     &                     *Coupl23_A* Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))) )

!(p2(Q)-radiation)*(p2(H)-radiation)

	 treeQS =  treeQS - 2*( (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc
     &                     *V2rad)*Coupl23_A
     &                     *Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))   )


!**********************************************************Dipoles( q*s(Higgs) )***************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeQS = treQS - 

!p2-splitting

!        treeQS = treQS -  

!     endif
  
!**********************************************************************************************************************************

         enddo

!Vector:- from Polarisationvector -(Propagator)*-(i-combination)=-   !Integrationsfehler riesig für vegas (~10) aber ok für Cuhre (~10^-3)

         do j= sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) 
            Coupl23_Bc = dconjg(JLR(j,argyh,argxh))
            propac = dconjg(DenSmV23(j)) 

!(p1(Q)-radiation)*(p1(V)-radiation)

	 treeQS =  treeQS - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))
     &                     *V1radc*V1rad)*Coupl23_A
     &                     *Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(Q)-radiation)*(p2(V)-radiation)+(p2(V)-radiation)*(p1(Q)-radiation)

	 treeQS =  treeQS - 2*( (sqrtS**2*(-2*p1k1 + 2*p1k2 + (muu3**2 - muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc
     &                     *V1rad)*Coupl23_A*Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(Q)-radiation)*(p1(V)-radiation)+(p1(V)-radiation)*(p2(Q)-radiation)

	 treeQS =  treeQS - 2*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                      V2rad)*Coupl23_A*Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   ) 

!(p2(Q)-radiation)*(p2(V)-radiation)

	 treeQS =  treeQS + 2*( (((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))
     &                     *V2radc*V2rad)*Coupl23_A
     &                     *Coupl23_Ac*Coupl23_Bc/
     &                      (propac*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))     )

!**********************************************************Dipoles( q*s(Vector) )**************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeQS = treeQS - 

!p2-splitting

!        treeQS = treeQS -  

!     endif
  
!**********************************************************************************************************************************

         enddo
                  
      QQ2xx_M2Tree23HH= QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQS)

!    	print *, "Tree_Mqs=", treeQS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!==============================================================QbQb================================================================
!- from Polarisationvector
    
      do i=sumstartV, sumendV 

         treeQbQb = (0d0,0d0)
         treeQbQ = (0d0,0d0)
         treeQbTU = (0d0,0d0)
         treeQbS = (0d0,0d0)

         Coupl23_A = VVSfSf(ztype,i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)
         Coupl23_B = JLR(i,argyh,argxh)
         propa = DenSmV23(i)    
  

         do j=sumstartV, sumendV 

            Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))
            Coupl23_Bc = dconjg(JLR(j,argyh,argxh)) 
            propac = dconjg(DenSmV23(j))     
 
         treeQbQb = treeQbQb -((sqrtS**2*(1 + 2*muu3**2 + 2*muu4**2 - muu5**2 - x1 - x2)*Coupl23_Bc*Coupl23_Ac*
     &                       Coupl23_B*Coupl23_A)/(propac*propa)) 
         
         enddo  


         QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQbQb)

!        print *, "Tree_Mqbqb=", treeQbQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!==============================================================QbQ=================================================================
!- from Polarisationvector - from Vector Propagator

         Coupl23_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))  

!Qb-Q(p1-radiation)
 
         treeQbQ = treeQbQ + 2*(((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*
     &                      Coupl23_Ac*
     &                      V1radc*Coupl23_B*
     &                      Coupl23_A)/
     &                      (2.*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))) 

!Qb-Q(p2-radiation)
 
         treeQbQ = treeQbQ + 2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu3**2 + 2*muu4**2 - x1 + x2))*
     &                      Coupl23_Ac*
     &                      V2radc*Coupl23_B*
     &                      Coupl23_A)/
     &                      (2.*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))
         

         QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQbQ)

!        print *, "Tree_Mqbq=", treeQbQ
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!==============================================================QbS(Higgs,Vector)===================================================
          
!**************************************************************QbS(Higgs)********************************************************** 
!- from Polarisationvector - from Vector Propagator - from i-combination  

         do j=sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))  
            Coupl23_Bc = dconjg(ILR(j,argxh,argyh))
            propac = dconjg(DenSmH23(j)) 
          
!Qb-S(p1-radiation)
 
         treeQbS = treeQbS +  2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*Coupl23_Ac*Coupl23_Bc*
     &                       V1radc*Coupl23_B*Coupl23_A)/
     &                       (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )
     
!Qb-S(p2-radiation)
 
         treeQbS = treeQbS +  2*(((4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu3**2 + 2*muu4**2 - x1 + x2))*Coupl23_Ac*Coupl23_Bc*
     &                       V2radc*Coupl23_B*Coupl23_A)/
     &                       (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))) 

         enddo

!**************************************************************QbS(Vector)*********************************************************
!- from Polarisationvector (- from Vector Propagator)^2 - from i-combination 

         do j=sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) 
            Coupl23_Bc = dconjg(JLR(j,argyh,argxh))
            propac = dconjg(DenSmV23(j)) 
          
!Qb-S(p1-radiation)
 
         treeQbS = treeQbS + 2*(((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*Coupl23_Bc*
     &                      Coupl23_Ac*V1radc*Coupl23_B*
     &                      Coupl23_A)/
     &                      (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!Qb-S(p2-radiation)
 
         treeQbS = treeQbS +  2*(-((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu3**2 + 2*muu4**2 - x1 + x2))*Coupl23_Bc*
     &                       Coupl23_Ac*V2radc*Coupl23_B*
     &                       Coupl23_A)/
     &                       (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

         enddo

!**********************************************************************************************************************************

         QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQbS)

!        print *, "Tree_Mqbs=", treeQbS
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!==================================================================QbT/U===========================================================
 
         do j= 1, 2

            Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
            Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
            sfmtc = MSf(j,proptypetc,propgentc)
            sfwtc = WSf(j,proptypetc,propgentc)

            Coupl23_Auc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
            Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
            sfmuc = MSf(j,proptypeuc,propgenuc)
            sfwuc = WSf(j,proptypeuc,propgenuc)     

            V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
            V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))          

!******************************************************************QbT*************************************************************         
!- from Polarisationvector - from Vector Propagator - from i-combination 
          
!Qb-T(p1-radiation)
 
         treeQbTU = treeQbTU + 2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*Coupl23_Atc*
     &                        Coupl23_Btc*V1radc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )))
     
!Qb-T(p2-radiation)
 
         treeQbTU = treeQbTU + 2*(((4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu3**2 + 2*muu4**2 - x1 + x2))*Coupl23_Atc*
     &                        Coupl23_Btc*V2radc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m1**2 + m3**2 - 2*p1k1 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

!Qb-T(p3-radiation)
 
         treeQbTU = treeQbTU + 2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(-2 - 4*muu3**2 - 4*muu4**2 + 2*muu5**2 + x1 + 3*x2))*Coupl23_Atc*
     &                        Coupl23_Btc*V3tradc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m1**2 + m3**2 - 2*p1k1 - sfmtc**2-cI*sfmtc*sfwtc)*(m2**2 + m4**2 + 2*p1k2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc))) 

!*******************************************************************QbU************************************************************ 
!- from Polarisationvector - from Vector Propagator - from i-combination 
         
!Qb-U(p1-radiation)
 
         treeQbTU = treeQbTU + 2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*Coupl23_Buc*
     &                        Coupl23_Auc*V1radc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )))

!Qb-U(p2-radiation)
 
         treeQbTU = treeQbTU + 2*(((4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu3**2 + 2*muu4**2 - x1 + x2))*Coupl23_Buc*
     &                        Coupl23_Auc*V2radc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m1**2 + m4**2 - 2*p1k2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

!Qb-U(p3-radiation)
 
         treeQbTU = treeQbTU + 2*(-((4*p1k1 - 4*p1k2 + sqrtS**2*(2 + 4*muu3**2 + 4*muu4**2 - 2*muu5**2 - 3*x1 - x2))*Coupl23_Buc*
     &                        Coupl23_Auc*V3uradc*Coupl23_B*
     &                        Coupl23_A)/
     &                        (2.*propa*(m1**2 + m4**2 - 2*p1k2 - sfmuc**2-cI*sfmuc*sfwuc)*(m2**2 + m3**2 + 2*p1k1 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)))

!**********************************************************************************************************************************

         enddo  

         QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeQbTU)

!        print *, "Tree_Mqbtu=", colorfactor23*dreal(treeQbTU)
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!**********************************************************************************************************************************

      enddo

!=====================================================================SS===========================================================

!************************************************************Higgs-Higgs (s-channel)***********************************************
!- from Polarisationvector

      treeSS = (0d0,0d0)

      do i= sumstartH, sumendH      

         Coupl23_A = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
         Coupl23_B = ILR(i,argxh,argyh)
         propa = DenSmH23(i)
                            
         do j= sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)) * onshell_H(j) 
            Coupl23_Bc = dconjg(ILR(j,argxh,argyh))
            propac = dconjg(DenSmH23(j)) 


          
!(p1-radiation)*(p1-radiation)

	 treeSS = treeSS - ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))

!(p1-radiation)*(p2-radiation)+(p2-radiation)*(p1-radiation)


	 treeSS = treeSS + 2d0*(   ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2
     &                     )))        )

!(p2-radiation)*(p2-radiation)

	 treeSS = treeSS + (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_A*Coupl23_Ac
     &                   *Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - sqrtS**2*x1 - sqrtS**2*x2 - m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

!**********************************************************Dipoles( s(Higgs)*s(Higgs) )********************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treSS - 

!p2-splitting

!        treeSS = treSS -  

!     endif
 
!**********************************************************************************************************************************
           
         enddo  
      enddo

!****************************************************Vector-Vector (s-Kanal)*******************************************************
!- from Polarisationvector
    
      do i=sumstartV, sumendV 

         Coupl23_A = VSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
         Coupl23_B = JLR(i,argyh,argxh)
         propa = DenSmV23(i)

         do j=sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))  * onshell_V(j) 
            Coupl23_Bc = dconjg(JLR(j,argyh,argxh))
            propac = dconjg(DenSmV23(j)) 
          
!(p1-radiation)*(p1-radiation)

         treeSS = treeSS - ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))**2
     &                   *V1radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))
    
!(p1-radiation)*(p2-radiation)+(p2-radiation)*(p1-radiation)

         treeSS = treeSS - 2d0*( (sqrtS**2*(-2*p1k1 + 2*p1k2 + (muu3**2 - muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*
     &                    V2radc*V1rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))    )


!(p2-radiation)*(p2-radiation)

         treeSS = treeSS + (((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)**2*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                    V2radc*V2rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))) 


!**********************************************************Dipoles( s(Vector)*s(Vector) )******************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treSS - 

!p2-splitting

!        treeSS = treSS -  

!     endif
 
!**********************************************************************************************************************************

         enddo


!**********************************************************Vector-Higgs (s-channel)************************************************
!- from Vector-Propagator *  - from Polarisationvector 

         do j=sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))  
            Coupl23_Bc = dconjg(ILR(j,argxh,argyh))
            propac = dconjg(DenSmH23(j)) 
          
!(p1-radiation)*(p1-radiation) 

         treeSS = treeSS + 2d0*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))   
     &                   *V1radc*
     &                    V1rad)* Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))     )

!(p1(V)-radiation)*(p2(H)-radiation)+(p2(H)-radiation)*(p1(V)-radiation)

         treeSS = treeSS + 2d0*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V2radc*
     &                    V1rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))    )

!(p2(V)-radiation)*(p1(H)-radiation)+(p1(H)-radiation)*(p2(V)-radiation)

         treeSS = treeSS + 2d0*( (sqrtS**2*(-2*p1k1 + 2*p1k2 + (muu3**2 - muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                    V2rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc /
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p2-radiation)*(p2-radiation) 

         treeSS = treeSS - 2d0*(  (((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))
     &                   *V2radc*V2rad)*Coupl23_A*Coupl23_Ac
     &                   * Coupl23_B*Coupl23_Bc/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))    )

!**********************************************************Dipoles( s(Vector)*s(Higgs) + s(Higgs)*s(Vector))***********************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS = treSS - 

!p2-splitting

!        treeSS = treSS -  

!     endif
 
!**********************************************************************************************************************************

         enddo     
      enddo

!**********************************************************************************************************************************

      QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeSS)

!        print *, "Tree_Mss=", treeSS
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH

!======================================================sFermion (t- and u-channel)=================================================

!-from Polarisationvector

      do i=1, 2  

         treeTU = (0d0,0d0) 
         treeTUQ = (0d0,0d0)
         treeTUS = (0d0,0d0) 

         Coupl23_At = HSfSft(argxh,sfmeigen1,i,sftype1,proptypet,sfgen1,propgent)
         Coupl23_Bt = HSfSfu(argyh,i,sfmeigen2,proptypet,sftype2,propgent,sfgen2)
         sfmt = MSf(i,proptypet,propgent)
         sfwt = WSf(i,proptypet,propgent)

         Coupl23_Au = HSfSfu(argyh,sfmeigen1,i,sftype1,proptypeu,sfgen1,propgenu)
         Coupl23_Bu = HSfSft(argxh,i,sfmeigen2,proptypeu,sftype2,propgenu,sfgen2)
         sfmu = MSf(i,proptypeu,propgenu)
         sfwu = WSf(i,proptypeu,propgenu)

         V3trad = VSfSf(ztype,i,i,proptypet,proptypet,propgent,propgent)                      !t-channel-radiation
         V3urad = VSfSf(ztype,i,i,proptypeu,proptypeu,propgenu,propgenu)                      !u-channel-radiation   
                      
         do j=1, 2         

           Coupl23_Atc = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
           Coupl23_Btc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
           sfmtc = MSf(j,proptypetc,propgentc)
           sfwtc = WSf(j,proptypetc,propgentc)

           Coupl23_Auc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
           Coupl23_Buc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
           sfmuc = MSf(j,proptypeuc,propgenuc)
           sfwuc = WSf(j,proptypeuc,propgenuc)     

           V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
           V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))             

!********************************************************(t*t)*********************************************************************
!(p1-radiation)*(p1-radiation) 

         treeTU = treeTU - ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))


!(p2-radiation)*(p2-radiation) 

         treeTU = treeTU + (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 - 2*p1k1 + m1**2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))  

!(channel-radiation)*(channel-radiation) 

         treeTU = treeTU + (((4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu1**2 - 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + 2*x2))*
     &                    V3tradc*V3trad)*Coupl23_At
     &                   *Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 - 2*p1k1 + m1**2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc)))

!(p1-radiation)*(p2-radiation)+(p2-radiation)*(p1-radiation)
 
         treeTU = treeTU + 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 - 2*p1k1 + m1**2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))         )

!(p1-radiation)*(channel-radiation)+(channel-radiation)*(p1-radiation)
 
         treeTU = treeTU - 2*( ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*V3tradc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 - 2*p1k1 + m1**2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!(p2-radiation)*(channel-radiation)+(channel-radiation)*(p2-radiation)
 
         treeTU = treeTU - 2*( ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V3tradc*
     &                    V2rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Atc*Coupl23_Btc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m3**2 - 2*p1k1 + m1**2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - sfmtc**2-cI*sfmtc*sfwtc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))      )

!*******************************************************Dipoles (t*t)**************************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU -  

!     endif
 
!*********************************************************(t*u+u*t)***************************************************************

!(p1(t)-radiation)*(p1(u)-radiation)+(p1(u)-radiation)*(p1(t)-radiation)

         treeTU = treeTU - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(p2(t)-radiation)*(p2(u)-radiation)+(p2(u)-radiation)*(p2(t)-radiation)

         treeTU = treeTU + 2*( (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))      )

!(channel(t)-radiation)*(channel(u)-radiation)+(channel(u)-radiation)*(channel(t)-radiation)

         treeTU = treeTU - 2*(  (sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2)*V3uradc*
     &                    V3trad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc))      )

!(p1(t)-radiation)*(p2(u)-radiation)+(p2(u)-radiation)*(p1(t)-radiation)

         treeTU = treeTU + 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))      )

!(p2(t)-radiation)*(p1(u)-radiation)+(p1(u)-radiation)*(p2(t)-radiation)

         treeTU = treeTU + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                    V2rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))   )

!(p1(t)-radiation)*(channel(u)-radiation)+(channel(u)-radiation)*(p1(t)-radiation)

         treeTU = treeTU - 2*(  ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*V3uradc*
     &                    V1rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(channel(t)-radiation)*(p1(u)-radiation)+(p1(u)-radiation)*(channel(t)-radiation)

         treeTU = treeTU - 2*( ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*V1radc*
     &                    V3trad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p2(t)-radiation)*(channel(u)-radiation)+(channel(u)-radiation)*(p2(t)-radiation)

         treeTU = treeTU - 2*( ((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V3uradc*
     &                    V2rad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))      )

!(channel(t)-radiation)*(p2(u)-radiation)+(p2(u)-radiation)*(channel(t)-radiation)

         treeTU = treeTU - 2*( ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V2radc*
     &                    V3trad)*Coupl23_At*Coupl23_Bt
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))      )
    
!****************************************************Dipoles(t*u+u*t)**************************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU - 

!     endif
                                  
!**********************************************************(u*u)*******************************************************************

!(p1-radiation)*(p1-radiation)

         treeTU = treeTU - ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad)*Coupl23_Au*Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!(p2-radiation)*(p2-radiation)

         treeTU = treeTU + (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad)*Coupl23_Au*Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

!(channel-radiation)*(channel-radiation)

         treeTU = treeTU - ((4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu1**2 + 2*muu2**2 + 2*muu3**2 + 2*muu4**2 - muu5**2 - 2*x1))*
     &                    V3uradc*V3urad)*Coupl23_Au
     &                   *Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc))  

!(p1-radiation)*(p2-radiation)+(p2-radiation)*(p1-radiation)

         treeTU = treeTU + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad)*Coupl23_Au*Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                    sqrtS**2*x1 - sqrtS**2*x2  )))       )
  
!(p1-radiation)*(channel-radiation)+(channel-radiation)*(p1-radiation)

         treeTU = treeTU - 2*( ((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*V3uradc*
     &                    V1rad)*Coupl23_Au*Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                    (m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))       )

!(p2-radiation)*(channel-radiation)+(channel-radiation)*(p2-radiation)

         treeTU = treeTU - 2*( ((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V3uradc*
     &                    V2rad)*Coupl23_Au*Coupl23_Bu
     &                   *Coupl23_Auc*Coupl23_Buc
     &                   /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m4**2 - 2*p1k2 + m1**2 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmuc**2-cI*sfmuc*sfwuc)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!*******************************************************Dipoles(u*u)***************************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTU = treeTU - 

!p2-splitting

!        treeTU = treeTU - 

!     endif
 
!**********************************************************************************************************************************

         enddo

         QQ2xx_M2Tree23HH = QQ2xx_M2Tree23HH + colorfactor23*dreal(treeTU)

!    	print *, "Tree_Mtt+Tree_Muu+Tree_Mut=", treeTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH 

!====================================================t*q + u*q (sFermion)==========================================================
!- from Polarisationvector * - from i-combination    

         Coupl23_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))        

!********************************************************( t*q )*******************************************************************
!(p1(t)-radiation)*(p1(q)-radiation)

         treeTUQ =  treeTUQ + 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                        V1rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p1(t)-radiation)*(p2(q)-radiation)

         treeTUQ =  treeTUQ - 2*(  ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                        V1rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  )))                 )

!(p2(t)-radiation)*(p1(q)-radiation)

         treeTUQ =  treeTUQ - 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                        V2rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  )))               )

!(p2(t)-radiation)*(p2(q)-radiation)

         treeTUQ =  treeTUQ - 2*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))                  )

!(channel(t)-radiation)*(p1(q)-radiation)

         treeTUQ =  treeTUQ + 2*( ((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*V1radc*
     &                        V3trad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))                )

!(channel(t)-radiation)*(p2(q)-radiation)

         treeTUQ =  treeTUQ + 2*(  ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V2radc*
     &                        V3trad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac
     &                       /((m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))                 )

!******************************************************Dipoles( t*q )**************************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUQ = treeTUQ - 

!p2-splitting

!        treeTUQ = treeTUQ -  

!     endif
  
!**********************************************************( u*q )*****************************************************************
!(p1(u)-radiation)*(p1(q)-radiation)  

         treeTUQ =  treeTUQ + 2*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc
     &                       *V1rad)*Coupl23_Au*Coupl23_Bu
     &                       *Coupl23_Ac
     &                       /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))        )     

!(p1(u)-radiation)*(p2(q)-radiation) 

        treeTUQ =  treeTUQ - 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad)
     &                      *Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac
     &                      /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))       )     

!(p2(u)-radiation)*(p1(q)-radiation)  

         treeTUQ =  treeTUQ - 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*V2rad)
     &                       *Coupl23_Au*Coupl23_Bu*Coupl23_Ac
     &                       /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))       )     

!(p2(u)-radiation)*(p2(q)-radiation)  

         treeTUQ =  treeTUQ - 2*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)
     &                       *Coupl23_Au*Coupl23_Bu*Coupl23_Ac
     &                       /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))   )     

!(channel(u)-radiation)*(p1(q)-radiation)  

         treeTUQ =  treeTUQ + 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*V1radc*V3urad)
     &                       *Coupl23_Au*Coupl23_Bu                                                                                                                       
     &                       *Coupl23_Ac
     &                       /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))       )     

!(channel(u)-radiation)*(p2(q)-radiation)  

         treeTUQ =  treeTUQ + 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V2radc*V3urad)
     &                       *Coupl23_Au*Coupl23_Bu
     &                       *Coupl23_Ac
     &                       /((m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )     

!******************************************************Dipoles( u*q )**************************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUQ = treeTUQ - 

!p2-splitting

!        treeTUQ = treeTUQ -  

!     endif
  
!**********************************************************************************************************************************            

         QQ2xx_M2Tree23HH= QQ2xx_M2Tree23HH + colorfactor23*dreal(treeTUQ)

!    	print *, "Tree_Mtq+Tree_Muq=", treeTUQ
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23HH
    

!=============================================u*s (Higgs)+ u*s (Vector) + t*s (Higgs)+ t*s (Vector)================================

!***************************************************************( t*s(Higgs) )*****************************************************

         do j=sumstartH, sumendH 

            Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))  
            Coupl23_Bc = dconjg(ILR(j,argxh,argyh))
            propac = dconjg(DenSmH23(j)) 

!- from Polarisationvector

!(p1(t)-radiation)*(p1(s)-radiation)

         treeTUS =  treeTUS - 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc
     &                    *V1rad)*Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))      )

!(p1(t)-radiation)*(p2(s)-radiation)

         treeTUS =  treeTUS + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad)
     &                    *Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))     )

!(p2(t)-radiation)*(p1(s)-radiation)

         treeTUS =  treeTUS + 2*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*V2rad)
     &                    *Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))      )

!(p2(t)-radiation)*(p2(s)-radiation)

         treeTUS =  treeTUS + 2*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)
     &                    *Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))      )

!(p3(t)-radiation)*(p1(s)-radiation)

         treeTUS =  treeTUS - 2*(((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*V1radc*V3trad)
     &                    *Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))     )

!(p3(t)-radiation)*(p2(s)-radiation)

         treeTUS =  treeTUS - 2*( ((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V2radc*V3trad)
     &                    *Coupl23_At*Coupl23_Bt
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))     )

!******************************************************Dipoles( t*s(Higgs) )*******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUS = treeTUS - 

!p2-splitting

!        treeTUS = treeTUS -  

!     endif
  
!**********************************************************( u*s(Higgs) )**********************************************************    
!(p1(u)-radiation)*(p1(s)-radiation)         

         treeTUS =  treeTUS - 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*V1rad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))    )         

!(p1(u)-radiation)*(p2(s)-radiation)


         treeTUS =  treeTUS + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))    )         

!(p2(u)-radiation)*(p1(s)-radiation)

         treeTUS =  treeTUS + 2*( ((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*V2rad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))    )         

!(p2(u)-radiation)*(p2(s)-radiation)

         treeTUS =  treeTUS + 2*( (((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))   )
         
!(p3(u)-radiation)*(p1(s)-radiation)

         treeTUS =  treeTUS - 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*V1radc*V3urad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))    )         

!(p3(u)-radiation)*(p2(s)-radiation)
    
         treeTUS =  treeTUS - 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V2radc*V3urad)
     &                    *Coupl23_Au*Coupl23_Bu
     &                    *Coupl23_Ac*Coupl23_Bc
     &                    /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                     (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )         

!******************************************************Dipoles( t*s(Higgs) )*******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUS = treeTUS - 

!p2-splitting

!        treeTUS = treeTUS -  

!     endif
  
!**********************************************************( t*s(Vector) )*********************************************************    
         
         enddo
         
         do j=sumstartV, sumendV 

            Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) 
            Coupl23_Bc = dconjg(JLR(j,argyh,argxh))
            propac = dconjg(DenSmV23(j)) 

! - form Propagator * - from Polarisationvector 

!(p1(t)-radiation)*(p1(s)-radiation)

	 treeTUS =  treeTUS + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                        V1rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p1(t)-radiation)*(p2(s)-radiation)

	 treeTUS =  treeTUS + 2*((sqrtS**2*(-2*p1k1 + 2*p1k2 + (muu3**2 - muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                        V1rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p1(s)-radiation)

	 treeTUS =  treeTUS + 2*((sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                        V2rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p2(s)-radiation)

	 treeTUS =  treeTUS - 2*((((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p3(t)-radiation)*(p1(s)-radiation)

	 treeTUS =  treeTUS - 2*((((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*(-2*p1k1 + 2*p1k2 + sqrtS**2*(-muu3**2 + muu4**2 + x1 - x2))*
     &                        V1radc*V3trad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )))  )

!(p3(t)-radiation)*(p2(s)-radiation)

	 treeTUS =  treeTUS + 2*(((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*V2radc*
     &                        V3trad)*Coupl23_At*Coupl23_Bt
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(propac*(m3**2 - 2*p1k1 + m1**2 -sfmt**2+cI*sfmt*sfwt)*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 -sfmt**2+cI*sfmt*sfwt)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!******************************************************Dipoles( t*s(Vectors) )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUS = treeTUS - 

!p2-splitting

!        treeTUS = treeTUS -  

!     endif
  
!**********************************************************( u*s(Vector) )*********************************************************    

!(p1(u)-radiation)*(p1(s)-radiation)

        treeTUS =  treeTUS + 2*( ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                       V1rad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

!(p1(u)-radiation)*(p2(s)-radiation)

        treeTUS =  treeTUS + 2*( (sqrtS**2*(-2*p1k1 + 2*p1k2 + (muu3**2 - muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                       V1rad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p2(u)-radiation)*(p1(s)-radiation)

        treeTUS =  treeTUS + 2*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                       V2rad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p2(u)-radiation)*(p2(s)-radiation)

        treeTUS =  treeTUS - 2*( (((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                       V2rad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p3(u)-radiation)*(p1(s)-radiation)

        treeTUS =  treeTUS + 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*V1radc*
     &                       V3urad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!(p3(u)-radiation)*(p2(s)-radiation)

        treeTUS =  treeTUS + 2*(((2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*V2radc*
     &                       V3urad)*Coupl23_Au*Coupl23_Bu
     &                      *Coupl23_Ac*Coupl23_Bc
     &                      /(propac*(m4**2 - 2*p1k2 + m1**2 - sfmu**2+cI*sfmu*sfwu)*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - sfmu**2+cI*sfmu*sfwu)*
     &                       (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   )

!******************************************************Dipoles( u*s(Higgs) )*******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeTUS = treeTUS - 

!p2-splitting

!        treeTUS = treeTUS -  

!     endif
 
!**********************************************************************************************************************************    

         enddo

         QQ2xx_M2Tree23HH= QQ2xx_M2Tree23HH + colorfactor23*dreal(treeTUS)
 
!    	print *, "Tree_Mts+Tree_Mus=", treeTUS
!     	print *, "Tree_M2Tree23HH=", QQ2xx_M2Tree23HH

!==================================================================================================================================
          
      enddo

!       QQ2xx_M2Tree23HH = treeTU   !Nur bestimmte Kanaele sollen zaehlen

       end
