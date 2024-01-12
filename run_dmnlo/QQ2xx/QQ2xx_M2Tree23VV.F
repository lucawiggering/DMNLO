!**********************************************************************************************************************************

	double precision function QQ2xx_M2Tree23VV()

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to3.h"

	

	integer i,j
        double precision sfmt,sfmu,sfmtc,sfmuc,sfwt,sfwu,sfwtc,sfwuc
	complex*16 tree23QQ,tree23QbQb,tree23QcQc,tree23QbQ,tree23QcQ,tree23QcQb,tree23SS,
     &             tree23QS,tree23QbS,tree23QcS,tree23TU,tree23QTU,tree23QbTU,tree23QcTU,tree23STU,
     &             Coupl23_A,Coupl23_At,Coupl23_Au,Coupl23_Ac,Coupl23_Atc,Coupl23_Auc,Coupl23_B,Coupl23_Bt,Coupl23_Bu,
     &             Coupl23_Bc,Coupl23_Btc,Coupl23_Buc,Coupl23_Ct,Coupl23_Cu,Coupl23_Ctc,Coupl23_Cuc,Coupl23_Dt,
     &             Coupl23_Du,Coupl23_Dtc,Coupl23_Duc,propa,propac
        
        QQ2xx_M2Tree23VV = 0d0    

! =================================================================================================================================	
!                         ---               Vector-Vector-final-states (+ Ghost-Ghost-final-states)                   ---
! =================================================================================================================================

! --- Mqq ---

!===============================================================q*q================================================================
    
!Factor 4 from polarisation-vector-contraction * -from Vector-contraction

        tree23QQ = 0d0
        tree23QS = 0d0

        Coupl23_A = VVSfSf(xtype,ytype,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2)
        Coupl23_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))

!(p1-radiation)*(p1-radiation)
 
        tree23QQ = tree23QQ - 4*((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                        V1rad)*Coupl23_A*Coupl23_Ac/
     &                        ((2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!(p1-radiation)*(p2-radiation)
!additional (-1) from amplitude
 
        tree23QQ = tree23QQ + 8*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                        V1rad)*Coupl23_A*Coupl23_Ac/
     &                        ((2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  )))    )

!(p2-radiation)*(p2-radiation)
!additional (-1) from amplitude
 
        tree23QQ = tree23QQ + 4*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_A*Coupl23_Ac/
     &                        ((2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 - 
     &                        m2**2 + m2**2)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))


!**********************************************************Dipoles( q*q )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree223QQ = tree23QQ - 

!p2-splitting

!        tree23QQ = tree23QQ -  

!      endif
  
!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QQ)

!        print *, "Tree_Mqq=", colorfactor23*tree23QQ
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Msq+Mqs ---

!===========================================================q*s (Higgs)+ q*s (Vector)==============================================

! - from i-combination  * -from Vector-contraction (=0 for Gluon-final-state wire couplings)
     
        do j=sumstartH, sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_Bc = dconjg(LLR(xtype,ytype,j))
          propac = dconjg(DenSmH23(j))

!(p1(q)-radiation)*(p1(Higgs)-radiation)

	tree23QS =  tree23QS + 8*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                          V1rad)*Coupl23_Ac
     &                         *Coupl23_Bc*Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )) )

!(p2(q)-radiation)*(p1(Higgs)-radiation)


	tree23QS =  tree23QS - 8*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                          V2rad)*Coupl23_Ac
     &                         *Coupl23_Bc*Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                          sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p1(q)-radiation)*(p2(Higgs)-radiation)
! (-1) from amplitude

	tree23QS =  tree23QS - 8*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                          V1rad)*Coupl23_Ac
     &                         *Coupl23_Bc*Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                          sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p2(q)-radiation)*(p2(Higgs)-radiation)


	tree23QS =  tree23QS - 8*((((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                          V2rad)*Coupl23_Ac
     &                         *Coupl23_Bc*Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                          sqrtS**2*x2  )*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))  )

!**********************************************************Dipoles( Higgs*q )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QS = tree23QS - 

!p2-splitting

!        tree23QS = tree23QS -  

!      endif
  
!**********************************************************************************************************************************

        enddo

!-(Propagator)*-(i-combination) * (-1) from Vector-contraction   

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))

!(p1(q)-radiation)*(p1(Vector)-radiation)

	tree23QS =  tree23QS - 2*((-3*(4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*
     &                          V1radc*V1rad)
     &                         *Coupl23_Ac*Coupl23_Bc
     &                         *Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )) )

!(p2(q)-radiation)*(p1(Vector)-radiation)

	tree23QS =  tree23QS - 2*((-3*sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*
     &                          V1radc*V2rad)
     &                         *Coupl23_Ac*Coupl23_Bc
     &                         *Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                          sqrtS**2*x1 - sqrtS**2*x2  )) )


!(p1(q)-radiation)*(p2(Vector)-radiation)

	tree23QS =  tree23QS - 2*((3*sqrtS**2*(2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                          V2radc*V1rad)
     &                         *Coupl23_Ac*Coupl23_Bc
     &                         *Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                          sqrtS**2*x1 - sqrtS**2*x2  )) )

!(p2(q)-radiation)*(p2(Vector)-radiation)

	tree23QS =  tree23QS - 2*((3*(2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*
     &                          (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                          V2rad)*Coupl23_Ac
     &                         *Coupl23_Bc*Coupl23_A/
     &                          (propac*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )) )



!**********************************************************Dipoles( Vector*q )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QS = tree23QS - 

!p2-splitting

!        tree23QS = tree23QS -  

!      endif
  
!**********************************************************************************************************************************

        enddo

        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QS)

!    	print *, "Tree_Mqs=", colorfactor23*tree23QS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mqbqb ---
! qb means additional quartic coupling through gluon radiation i.e. tT->Z->WW will have tTGZ coupling in 2->3
!==========================================================================qb*qb===================================================
   
      do i=sumstartV, sumendV
     
        tree23QbQb = 0d0
        tree23QbQ = 0d0
        tree23QbS = 0d0
        tree23QbTU = 0d0

        Coupl23_A = VVSfSf(ztype,i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2)
        Coupl23_B = VVV(i,xtype,ytype)
        propa = DenSmV23(i)

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2)) * onshell_V(j)
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))
 
! -from Vector-contraction 
 
        tree23QbQb = tree23QbQb + (-9*sqrtS**2*(-1 + muu3**2 + muu4**2 + muu5**2 + x1 + x2)*Coupl23_Ac*Coupl23_Bc*
     &                             Coupl23_A*Coupl23_B)/(propac*propa)

!Ghosts

        tree23QbQb = tree23QbQb + (sqrtS**2*(1 + muu3**2 + muu4**2 - muu5**2 - x1 - x2)*dconjg(GGV(xtype,ytype,j))*Coupl23_Ac*GGV(xtype,ytype,i)*
     &                            Coupl23_A)/(propac*propa)

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QbQb)

!        print *, "Tree_Mqbqb=", colorfactor23*tree23QbQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mqbq ---

!=======================================================================qb*q=======================================================
 
! -from Vector-contraction *  - from Propagator

        Coupl23_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))

!(p1-radiation)
 
        tree23QbQ = tree23QbQ - 3*(4*p1k1 - 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 - x1 + x2))*dreal(Coupl23_A*Coupl23_B*Coupl23_Ac*V1radc/
     &    (propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + muu2**2*sqrtS**2 + muu5**2*sqrtS**2 - m1**2)))

!(p2-radiation)

 
        tree23QbQ = tree23QbQ + 3*(-4*p1k1 + 4*p1k2 + sqrtS**2*(2*muu3**2 - 2*muu4**2 + x1 - x2))*dreal(-(Coupl23_A*Coupl23_B*Coupl23_Ac*V2radc/
     &      (propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - muu1**2*sqrtS**2 - muu5**2*sqrtS**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!**********************************************************Dipoles( qb*q )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree223QbQ = tree23QbQ - 

!p2-splitting

!        tree23QbQ = tree23QbQ -  

!      endif
  
!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QbQ)

!        print *, "Tree_Mqbq=", colorfactor23*tree23QbQ
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- MqbH+MqbV ---

!===================================================qb*s (Higgs)+ qb*s (Vector)====================================================
        

!************************************************************Higgs*****************************************************************

! - from i-combination  * -from Vector-contraction *- from Vectorpropagator

        do j=sumstartH, sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_Bc = dconjg(LLR(xtype,ytype,j))
          propac = dconjg(DenSmH23(j))

!(p1(Higgs)-radiation)

	tree23QbS =  tree23QbS - 2*(colorfactor23*(-3*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*Coupl23_Ac*Coupl23_Bc*
     &                            V1radc*Coupl23_A*Coupl23_B)/
     &                            (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!Ghosts

!(p1(Higgs)-radiation)

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

        tree23QbS = tree23QbS + 2*((colorfactor23*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*GGH(ytype,xtype,j)*Coupl23_Ac*
     &                           V1radc*GGV(xtype,ytype,i)*Coupl23_A)/
     &                           (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

        tree23QbS = tree23QbS +  2*(-(colorfactor23*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*
     &                            GGH(xtype,ytype,j)*HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)*V1rad)/
     &                            (2.*dconjg(propa)*DenSmH23(j)*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

        else

        tree23QbS = tree23QbS + 2*((colorfactor23*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*GGH(ytype,xtype,j)*Coupl23_Ac*
     &                           V1radc*GGV(xtype,ytype,i)*Coupl23_A)/
     &                           (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

        tree23QbS = tree23QbS +  2*(-(colorfactor23*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*
     &                            GGH(xtype,ytype,j)*HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)*V1rad)/
     &                            (2.*dconjg(propa)*DenSmH23(j)*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

        endif

!(p2(Higgs)-radiation)

	tree23QbS =  tree23QbS -2*(colorfactor23*(3*(2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*Coupl23_Ac*Coupl23_Bc*
     &                           V2radc*Coupl23_A*Coupl23_B)/
     &                           (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!Ghosts

!(p2(Higgs)-radiation)

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

        tree23QbS = tree23QbS + 2*(-(colorfactor23*(4*p1k2 - sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + 2*x2 + x1))*GGH(ytype,xtype,j)*Coupl23_Ac*
     &                           V2radc*GGV(xtype,ytype,i)*Coupl23_A)/
     &                           (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        tree23QbS = tree23QbS + 2*((colorfactor23*(4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 - x2 - 2*x1))*dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*
     &                           GGH(xtype,ytype,j)*HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)*V2rad)/
     &                           (2.*dconjg(propa)*DenSmH23(j)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        else

        tree23QbS = tree23QbS + 2*(-(colorfactor23*(4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*GGH(ytype,xtype,j)*Coupl23_Ac*
     &                           V2radc*GGV(xtype,ytype,i)*Coupl23_A)/
     &                           (2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        tree23QbS = tree23QbS + 2*((colorfactor23*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - 2*x2))*dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*
     &                           GGH(xtype,ytype,j)*HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)*V2rad)/
     &                           (2.*dconjg(propa)*DenSmH23(j)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        endif

!**********************************************************Dipoles( Higgs*qb )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QbS = tree23QbS - 

!p2-splitting

!        tree23QbS = tree23QbS -  

!      endif
  
!**********************************************************************************************************************************

        enddo

!********************************************************Vector********************************************************************
! -(i-combination) * -from Vector-contraction   

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j)
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))


!(p1(Vector)-radiation)

	tree23QbS =  tree23QbS + 2*(-(colorfactor23*(32*p1k1**2 + 32*p1k2**2 + 4*p1k2*sqrtS**2*(-2 + 6*muu3**2 - 14*muu4**2 + 2*muu5**2 - 19*x1 - 4*x2) + 
     &                           4*p1k1*(56*p1k2 - sqrtS**2*(2 + 14*muu3**2 - 6*muu4**2 - 2*muu5**2 + 4*x1 + 19*x2)) + 
     &                           sqrtS**4*(-14 + 8*muu3**2 - 10*muu3**4 + 8*muu4**2 + 20*muu3**2*muu4**2 - 10*muu4**4 + 12*muu5**2 - 4*muu3**2*muu5**2 - 4*muu4**2*muu5**2 + 2*muu5**4 + 15*x1 + 17*muu3**2*x1 - 
     &                           13*muu4**2*x1 + muu5**2*x1 + 3*x1**2 + 15*x2 - 13*muu3**2*x2 + 17*muu4**2*x2 + muu5**2*x2 + 26*x1*x2 + 3*x2**2 - 
     &                           4*muu1**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2)) - 12*muu2**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2))))*Coupl23_Ac*V1radc*Coupl23_Bc*Coupl23_A*Coupl23_B)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!Ghosts

!(p1(Vector)-radiation)

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

        tree23QbS = tree23QbS - 2*((colorfactor23*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*(-4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 + x1 - x2))*dconjg(GGV(xtype,ytype,j))*
     &                           Coupl23_Ac*V1radc*GGV(xtype,ytype,i)*
     &                           Coupl23_A)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

        tree23QbS = tree23QbS - 2*(-(colorfactor23*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*dconjg(GGV(xtype,ytype,i))*
     &                           dconjg(Coupl23_A)*onshell_V(j)*GGV(xtype,ytype,j)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)*
     &                           V1rad)/
     &                           (4.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

        else

        tree23QbS = tree23QbS - 2*((colorfactor23*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*(-4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 + x2))*dconjg(GGV(xtype,ytype,j))*
     &                           Coupl23_Ac*V1radc*GGV(xtype,ytype,i)*
     &                           Coupl23_A)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

        tree23QbS = tree23QbS - 2*(-(colorfactor23*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*dconjg(GGV(xtype,ytype,i))*
     &                           dconjg(Coupl23_A)*onshell_V(j)*GGV(xtype,ytype,j)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)*
     &                           V1rad)/
     &                           (4.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

        endif

!(p2(Vector)-radiation)

	tree23QbS =  tree23QbS + 2*(colorfactor23*(32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(56*p1k2 + sqrtS**2*(2 + 14*muu3**2 - 6*muu4**2 - 2*muu5**2 - 4*x1 - 9*x2)) - 
     &                           4*p1k2*sqrtS**2*(-2 + 6*muu3**2 - 14*muu4**2 + 2*muu5**2 + 9*x1 + 4*x2) - 
     &                           sqrtS**4*(14 - 8*muu3**2 + 10*muu3**4 - 8*muu4**2 - 20*muu3**2*muu4**2 + 10*muu4**4 - 12*muu5**2 + 4*muu3**2*muu5**2 + 4*muu4**2*muu5**2 - 2*muu5**4 - 11*x1 + 11*muu3**2*x1 + 
     &                           muu4**2*x1 - 5*muu5**2*x1 - 3*x1**2 - 11*x2 + muu3**2*x2 + 11*muu4**2*x2 - 5*muu5**2*x2 - 6*x1*x2 - 3*x2**2 + 12*muu1**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2)) + 
     &                           4*muu2**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2))))*Coupl23_Ac*V2radc*Coupl23_Bc*Coupl23_A*Coupl23_B)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 -m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))

!Ghosts

!(p2(Vector)-radiation)

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

        tree23QbS = tree23QbS - 2*((colorfactor23*(4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 - x1 - x2))*(4*p1k2 - sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*dconjg(GGV(xtype,ytype,j))*
     &                           Coupl23_Ac*V2radc*GGV(xtype,ytype,i)*
     &                           Coupl23_A)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        tree23QbS = tree23QbS - 2*(-(colorfactor23*(4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 - 2*x1 - x2))*(-4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*
     &                           dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*onshell_V(j)*GGV(xtype,ytype,j)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)*
     &                           V2rad)/
     &                           (4.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2 - 
     &                           m2**2)))

        else

        tree23QbS = tree23QbS - 2*((colorfactor23*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*(4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*dconjg(GGV(xtype,ytype,j))*
     &                           Coupl23_Ac*V2radc*GGV(xtype,ytype,i)*
     &                           Coupl23_A)/
     &                           (4.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

        tree23QbS = tree23QbS - 2*(-(colorfactor23*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - 2*x2))*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                           dconjg(GGV(xtype,ytype,i))*dconjg(Coupl23_A)*onshell_V(j)*GGV(xtype,ytype,j)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)*
     &                           V2rad)/
     &                           (4.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2 - 
     &                           m2**2)))

        endif

!**********************************************************Dipoles( Vector*qb )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QbS = tree23QbS - 

!p2-splitting

!        tree23QbS = tree23QbS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
         
        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + dreal(tree23QbS)

!    	print *, "Tree_Mqbs=", tree23QbS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mqbt+Mqbu ---

!====================================================qb*t + qb*u (sFermion)========================================================

        do j= 1,2        

          Coupl23_Atc = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl23_Dtc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmtc = MSf(j,proptypetc,propgentc)
          sfwtc = WSf(j,proptypetc,propgentc)

          Coupl23_Auc = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl23_Duc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmuc = MSf(j,proptypeuc,propgenuc)
          sfwuc = WSf(j,proptypeuc,propgenuc)

          V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
          V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))

! - from i-combination * -from Vector-contraction * - from Propagator

!(p1(t)-radiation)

        tree23QbTU =  tree23QbTU - 2*(((32*p1k1**2 + 32*p1k2**2 - 4*p1k2*sqrtS**2*(-2 + 2*muu1**2 - 2*muu2**2 + 8*muu3**2 + 6*x1 + 5*x2) + 
     &                              sqrtS**4*(-2 + 4*muu3**2 + 2*muu3**4 + 12*muu4**2 - 4*muu3**2*muu4**2 + 2*muu4**4 - 4*muu3**2*muu5**2 - 4*muu4**2*muu5**2 + 2*muu5**4 + 11*x1 + 3*muu3**2*x1 + muu4**2*x1 - 
     &                              3*muu5**2*x1 + x1**2 - 4*muu1**2*(-1 + 2*muu4**2 + muu5**2 + 2*x1) - 5*x2 + 9*muu3**2*x2 - 5*muu4**2*x2 + 5*muu5**2*x2 + 2*x1*x2 + 5*x2**2 - 
     &                              4*muu2**2*(-3 + 4*muu3**2 + 2*muu4**2 + 3*muu5**2 + 4*x1 + 2*x2)) + 4*p1k1*(8*p1k2 + sqrtS**2*(-6 + 2*muu1**2 - 2*muu2**2 - 4*muu3**2 - 4*muu4**2 + 4*muu5**2 - 3*x1 + 4*x2)))
     &                             *V1radc*Coupl23_Atc*
     &                              Coupl23_Dtc*Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                              (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(t)-radiation)

        tree23QbTU =  tree23QbTU - 2*(((32*p1k1**2 + 32*p1k2**2 + 4*p1k2*sqrtS**2*(6 + 2*muu1**2 - 2*muu2**2 + 4*muu3**2 + 4*muu4**2 - 4*muu5**2 - 8*x1 - 5*x2) + 
     &                              4*p1k1*(8*p1k2 + sqrtS**2*(-2 - 2*muu1**2 + 2*muu2**2 + 8*muu4**2 - 3*x1 + 2*x2)) - 
     &                              sqrtS**4*(2 - 12*muu3**2 - 2*muu3**4 - 4*muu4**2 + 4*muu3**2*muu4**2 - 2*muu4**4 + 4*muu3**2*muu5**2 + 4*muu4**2*muu5**2 - 2*muu5**4 + x1 + 5*muu3**2*x1 + 7*muu4**2*x1 - 
     &                              5*muu5**2*x1 - 3*x1**2 + x2 + 7*muu3**2*x2 + 5*muu4**2*x2 - 5*muu5**2*x2 - 6*x1*x2 - 3*x2**2 + 4*muu2**2*(-1 + 2*muu3**2 + muu5**2 + x1 + x2) + 
     &                              4*muu1**2*(-3 + 2*muu3**2 + 4*muu4**2 + 3*muu5**2 + x1 + 5*x2)))*V2radc*
     &                              Coupl23_Atc*Coupl23_Dtc*
     &                              Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (sfmtc**2+cI*sfmtc*sfwtc))*
     &                              (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!(channel(t)-radiation)

        tree23QbTU =  tree23QbTU - 2*(((32*p1k1**2 + 32*p1k2**2 - 4*p1k2*sqrtS**2*(-6 + 6*muu1**2 - 6*muu2**2 + 4*muu3**2 - 4*muu4**2 + 4*muu5**2 + 7*x2) - 
     &                              sqrtS**4*(-2 - 16*muu3**2 - 2*muu3**4 - 16*muu4**2 + 4*muu3**2*muu4**2 - 2*muu4**4 + 4*muu5**2 + 4*muu3**2*muu5**2 + 4*muu4**2*muu5**2 - 2*muu5**4 + 3*x1 - muu3**2*x1 + 
     &                              5*muu4**2*x1 + muu5**2*x1 + x1**2 + 4*muu2**2*(-2 + 3*muu3**2 + muu4**2 + 2*muu5**2 + x1) + 4*muu1**2*(-2 + muu3**2 + 3*muu4**2 + 2*muu5**2 - 3*x1 - 2*x2) + 15*x2 - 
     &                              3*muu3**2*x2 + 7*muu4**2*x2 - 7*muu5**2*x2 - 6*x1*x2 - 5*x2**2) + 4*p1k1*(8*p1k2 - sqrtS**2*(6 + 6*muu1**2 - 6*muu2**2 + 4*muu3**2 - 4*muu4**2 - 4*muu5**2 + x1 + 4*x2)))*
     &                              V3tradc*Coupl23_Atc*
     &                              Coupl23_Dtc*Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (sfmtc**2+cI*sfmtc*sfwtc))*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmtc**2+cI*sfmtc*sfwtc))))


!**********************************************************Dipoles( qb*t )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QbTU = tree23QbTU - 

!p2-splitting

!        tree23QbTU = tree23QbTU -  

!      endif
  
!**********************************************************************************************************************************
 
! - from i-combination * -from Vector-contraction * - from Propagator 

!(p1(u)-radiation)*(p1(q)-radiation)
     
        tree23QbTU =  tree23QbTU - 2*(-((32*p1k1**2 + 32*p1k2**2 + 4*p1k2*sqrtS**2*(-6 + 2*muu1**2 - 2*muu2**2 - 4*muu3**2 - 4*muu4**2 + 4*muu5**2 + 4*x1 - 3*x2) + 
     &                              sqrtS**4*(-2 + 12*muu3**2 + 2*muu3**4 + 4*muu4**2 - 4*muu3**2*muu4**2 + 2*muu4**4 - 4*muu3**2*muu5**2 - 4*muu4**2*muu5**2 + 2*muu5**4 - 5*x1 - 5*muu3**2*x1 + 9*muu4**2*x1 + 
     &                              5*muu5**2*x1 + 5*x1**2 + 11*x2 + muu3**2*x2 + 3*muu4**2*x2 - 3*muu5**2*x2 + 2*x1*x2 + x2**2 - 4*muu1**2*(-1 + 2*muu3**2 + muu5**2 + 2*x2) - 
     &                              4*muu2**2*(-3 + 2*muu3**2 + 4*muu4**2 + 3*muu5**2 + 2*x1 + 4*x2)) + 4*p1k1*(8*p1k2 - sqrtS**2*(-2 + 2*muu1**2 - 2*muu2**2 + 8*muu4**2 + 5*x1 + 6*x2)))*
     &                              V1radc*Coupl23_Duc*
     &                              Coupl23_Auc*Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                              (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))
 
!(p2(u)-radiation)*(p2(q)-radiation)
     
        tree23QbTU =  tree23QbTU - 2*(-((32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(8*p1k2 + sqrtS**2*(6 + 2*muu1**2 - 2*muu2**2 + 4*muu3**2 + 4*muu4**2 - 4*muu5**2 - 5*x1 - 8*x2)) - 
     &                              4*p1k2*sqrtS**2*(2 + 2*muu1**2 - 2*muu2**2 - 8*muu3**2 - 2*x1 + 3*x2) - 
     &                              sqrtS**4*(2 - 4*muu3**2 - 2*muu3**4 - 12*muu4**2 + 4*muu3**2*muu4**2 - 2*muu4**4 + 4*muu3**2*muu5**2 + 4*muu4**2*muu5**2 - 2*muu5**4 + x1 + 5*muu3**2*x1 + 7*muu4**2*x1 - 
     &                              5*muu5**2*x1 - 3*x1**2 + x2 + 7*muu3**2*x2 + 5*muu4**2*x2 - 5*muu5**2*x2 - 6*x1*x2 - 3*x2**2 + 4*muu2**2*(-1 + 2*muu4**2 + muu5**2 + x1 + x2) + 
     &                              4*muu1**2*(-3 + 4*muu3**2 + 2*muu4**2 + 3*muu5**2 + 5*x1 + x2)))*V2radc*
     &                              Coupl23_Duc*Coupl23_Auc*
     &                              Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (sfmuc**2+cI*sfmuc*sfwuc))*
     &                              (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!(channel(u)-radiation)*(p1(q)-radiation)
     
        tree23QbTU =  tree23QbTU - 2*(-((32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(8*p1k2 + sqrtS**2*(6 - 6*muu1**2 + 6*muu2**2 + 4*muu3**2 - 4*muu4**2 - 4*muu5**2 - 7*x1)) - 
     &                              4*p1k2*sqrtS**2*(6 + 6*muu1**2 - 6*muu2**2 - 4*muu3**2 + 4*muu4**2 - 4*muu5**2 + 4*x1 + x2) - 
     &                              sqrtS**4*(-2 - 16*muu3**2 - 2*muu3**4 - 16*muu4**2 + 4*muu3**2*muu4**2 - 2*muu4**4 + 4*muu5**2 + 4*muu3**2*muu5**2 + 4*muu4**2*muu5**2 - 2*muu5**4 + 15*x1 + 7*muu3**2*x1 - 
     &                              3*muu4**2*x1 - 7*muu5**2*x1 - 5*x1**2 + 4*muu1**2*(-2 + 3*muu3**2 + muu4**2 + 2*muu5**2 - 2*x1 - 3*x2) + 3*x2 + 5*muu3**2*x2 - muu4**2*x2 + muu5**2*x2 - 6*x1*x2 + x2**2 + 
     &                              4*muu2**2*(-2 + muu3**2 + 3*muu4**2 + 2*muu5**2 + x2)))*V3uradc*
     &                              Coupl23_Duc*Coupl23_Auc*
     &                              Coupl23_A*Coupl23_B)/
     &                              (4.*propa*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (sfmuc**2+cI*sfmuc*sfwuc))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc)))) 

!**********************************************************Dipoles( qb*u )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QbTU = tree23QbTU - 

!p2-splitting

!        tree23QbTU = tree23QbTU -  

!      endif
  
!**********************************************************************************************************************************
   
        enddo

        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QbTU)

!    	print *, "Tree_Mqbt+Tree_Mqbu=", colorfactor23*tree23QbTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

!==================================================================================================================================

      enddo

! --- Mqbt+Mqbu ---

!========================================================================qc*qc=====================================================

      do i=1,2

        tree23QcQc = 0d0
        tree23QcQ = 0d0
        tree23QcQb = 0d0
        tree23QcS = 0d0
        tree23QcTU = 0d0

        Coupl23_At = VSfSft(xtype,i,sfmeigen1,proptypet,sftype1,propgent,sfgen1)
        Coupl23_Bt = VVSfSf(ztype,ytype,sfmeigen2,i,sftype2,proptypet,sfgen2)
        Coupl23_Ct = VVSfSf(ztype,xtype,i,sfmeigen1,proptypet,sftype1,propgent)
        Coupl23_Dt = VSfSfu(ytype,sfmeigen2,i,sftype2,proptypet,sfgen2,propgent)
        sfmt = MSf(i,proptypet,propgent)
        sfwt = WSf(i,proptypet,propgent)

        Coupl23_Au = VSfSfu(ytype,i,sfmeigen1,proptypeu,sftype1,propgenu,sfgen1)
        Coupl23_Bu = VVSfSf(ztype,xtype,sfmeigen2,i,sftype2,proptypeu,sfgen2)
        Coupl23_Cu = VVSfSf(ztype,ytype,i,sfmeigen1,proptypeu,sftype1,propgenu)
        Coupl23_Du = VSfSft(xtype,sfmeigen2,i,sftype2,proptypeu,sfgen2,propgenu)
        sfmu = MSf(i,proptypeu,propgenu)
        sfwu = WSf(i,proptypeu,propgenu)

        do j=1,2

          Coupl23_Atc = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl23_Btc = dconjg(VVSfSf(ztype,ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2))
          Coupl23_Ctc = dconjg(VVSfSf(ztype,xtype,j,sfmeigen1,proptypetc,sftype1,propgentc))
          Coupl23_Dtc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmtc = MSf(j,proptypetc,propgentc)
          sfwtc = WSf(j,proptypetc,propgentc)

          Coupl23_Auc = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl23_Buc = dconjg(VVSfSf(ztype,xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2))
          Coupl23_Cuc = dconjg(VVSfSf(ztype,ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc))
          Coupl23_Duc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmuc = MSf(j,proptypeuc,propgenuc)
          sfwuc = WSf(j,proptypeuc,propgenuc)
 

! Factor 4 from polarisation-vector-contraction * -from Vector-contraction
! radiation from t-channel vertex V1 or V2 (additional quartic coupling in 2->3)
!*************************************************************************t*t******************************************************
!(V1-radiation)*(V1-radiation)
 
        tree23QcQc = tree23QcQc + (  (4*(4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*Coupl23_Atc*
     &                            Coupl23_Btc*Coupl23_At*
     &                            Coupl23_Bt)/
     &                            ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )

!(V1-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc - ((-4*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                            (-((Coupl23_Dtc*Coupl23_Ctc*
     &                            Coupl23_At*Coupl23_Bt)/
     &                            ((2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (sfmt**2-cI*sfmt*sfwt))*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmtc**2+cI*sfmtc*sfwtc))))))

!(V2-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc + (  (-4*(4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*Coupl23_Dtc*
     &                            Coupl23_Ctc*Coupl23_Dt*
     &                            Coupl23_Ct)/
     &                            ((2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmt**2-cI*sfmt*sfwt))*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmtc**2+cI*sfmtc*sfwtc))) ) 

!****************************************************t*u***************************************************************************
!(V1-radiation)*(V1-radiation)
 
        tree23QcQc = tree23QcQc - ((-4*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 8*muu1**2 - muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                            ( (Coupl23_Auc*Coupl23_Buc*
     &                            Coupl23_At*Coupl23_Bt)/
     &                            ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc)))))

!(V1-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc - 8*sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(-((Coupl23_Duc*
     &                            Coupl23_Cuc*Coupl23_At*
     &                            Coupl23_Bt)/
     &                            ((2*p1k1 - (muu1**2 + muu3**2)*sqrtS**2 + (sfmt**2-cI*sfmt*sfwt))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc)))))

!(V2-radiation)*(V1-radiation)
 
        tree23QcQc = tree23QcQc - 8*sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*(-((Coupl23_Auc*
     &                            Coupl23_Buc*Coupl23_Dt*
     &                            Coupl23_Ct)/
     &                            ((2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmt**2-cI*sfmt*sfwt))*(2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (sfmuc**2+cI*sfmuc*sfwuc)))))

!(V2-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc - ((4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 - muu4**2 + muu5**2 - x1 - x2))*
     &                            ((Coupl23_Duc*Coupl23_Cuc*
     &                            Coupl23_Dt*Coupl23_Ct)/
     &                            ((2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmt**2-cI*sfmt*sfwt))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc)))))

!****************************************************u*u***************************************************************************
!(V1-radiation)*(V1-radiation)
 
        tree23QcQc = tree23QcQc - 4*(  (-4*p1k2 + (4*muu1**2 + muu4**2)*sqrtS**2)*Coupl23_Auc*
     &                            Coupl23_Buc/
     &                            ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc)))*
     &                            Coupl23_Au*Coupl23_Bu  )

!(V1-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc - ((4*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - x1 + x2))*
     &                            (-((Coupl23_Duc*Coupl23_Cuc*
     &                            Coupl23_Au*Coupl23_Bu)/
     &                            ((2*p1k2 - (muu1**2 + muu4**2)*sqrtS**2 + (sfmu**2-cI*sfmu*sfwu))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc))))))

!(V2-radiation)*(V2-radiation)
 
        tree23QcQc = tree23QcQc + (  (-4*(4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*Coupl23_Duc*
     &                            Coupl23_Cuc*Coupl23_Du*
     &                            Coupl23_Cu)/
     &                            ((2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmu**2-cI*sfmu*sfwu))*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) - (sfmuc**2+cI*sfmuc*sfwuc)))  )

!**********************************************************************************************************************************

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QcQc)

!        print *, "Tree_Mqcqc=", colorfactor23*tree23QcQc
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

!===============================================================qc*q===============================================================
 
!-from Vector-contraction

        Coupl23_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))

!***************************************************************t*q****************************************************************
!(p1-radiation)*(p1-radiation)
 
        tree23QcQ = tree23QcQ - (4*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - muu5**2 - x2))*
     &                          ((V1radc*Coupl23_Ac*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

!(p2-radiation)*(p1-radiation)
 
        tree23QcQ = tree23QcQ - (4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*
     &                          ((V1radc*Coupl23_Ac*
     &                          Coupl23_Dt*Coupl23_Ct)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p1-radiation)*(p2-radiation)
 
        tree23QcQ = tree23QcQ - (-8*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &                          (-(V2radc*Coupl23_Ac*
     &                          Coupl23_At*Coupl23_Bt)/
     &                          ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!(p2-radiation)*(p2-radiation)
 
        tree23QcQ = tree23QcQ + ((4*p1k1 + sqrtS**2*(1 - 2*muu1**2 - 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1))*
     &                          (-((V2radc*Coupl23_Ac*
     &                          Coupl23_Dt*Coupl23_Ct)/
     &                          ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!*****************************************************u*q**************************************************************************
!(p1-radiation)*(p1-radiation)
 
        tree23QcQ = tree23QcQ - (4*p1k1 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - muu5**2 - x1))*
     &                          ((V1radc*Coupl23_Ac*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

!(p2-radiation)*(p1-radiation)
 
        tree23QcQ = tree23QcQ - (8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*
     &                          ((V1radc*Coupl23_Ac*
     &                          Coupl23_Du*Coupl23_Cu)/
     &                          ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

!(p1-radiation)*(p2-radiation)
 
        tree23QcQ = tree23QcQ - (-4*p1k1 - 8*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &                          (-(V2radc*Coupl23_Ac*
     &                          Coupl23_Au*Coupl23_Bu)/
     &                          ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!(p2-radiation)*(p2-radiation)
 
        tree23QcQ = tree23QcQ + ((4*p1k2 + sqrtS**2*(1 - 2*muu1**2 - 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - x2))*
     &                           (-((V2radc*Coupl23_Ac*
     &                           Coupl23_Du*Coupl23_Cu)/
     &                           ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!**********************************************************Dipoles( qc*q )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree223QcQ = tree23QcQ - 

!p2-splitting

!        tree23QcQ = tree23QcQ -  

!      endif

!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QcQ)

!        print *, "Tree_Mqcq=", colorfactor23*tree23QcQ
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mqbqb ---

!====================================================================qc*qb=========================================================

        do j=sumstartV, sumendV

          Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))
 
!- from Propagator * -from Vector-contraction

!*********************************************************************t*qb*********************************************************
!(p2-radiation)
 
        tree23QcQb = tree23QcQb - 6*(2*p1k1 + 4*p1k2 + sqrtS**2*(1 + muu4**2 - muu5**2 - x1 - x2))*(Coupl23_Ac*Coupl23_Bc*
     &                            Coupl23_At*Coupl23_Bt)/
     &                            (propac*(m3**2 - 2*p1k1 +m1**2 - (sfmt**2-cI*sfmt*sfwt)))

!(p1-radiation)
 
        tree23QcQb = tree23QcQb - 6*(-4*p1k1 - 2*p1k2 + sqrtS**2*(1 + muu3**2 - muu5**2 + x1))*(Coupl23_Ac*Coupl23_Bc*Coupl23_Dt*
     &                            Coupl23_Ct)/
     &                            (propac*(2*p1k2 + sqrtS**2*(muu2**2 + muu4**2 - x2) - (sfmt**2-cI*sfmt*sfwt)))


!**********************************************************************u*qb********************************************************
!(p2-radiation)
 
        tree23QcQb = tree23QcQb + 6*(4*p1k1 + 2*p1k2 + sqrtS**2*(1 + muu3**2 - muu5**2 - x1 - x2))*(Coupl23_Ac*Coupl23_Bc*Coupl23_Au*
     &                            Coupl23_Bu)/
     &                            (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu)))

!(p1-radiation)
 
        tree23QcQb = tree23QcQb + 6*(-2*p1k1 - 4*p1k2 + sqrtS**2*(1 + muu4**2 - muu5**2 + x2))*(Coupl23_Ac*Coupl23_Bc*
     &                            Coupl23_Du*Coupl23_Cu)/
     &                            (propac*(2*p1k1 + sqrtS**2*(muu2**2 + muu3**2 - x1) -(sfmu**2-cI*sfmu*sfwu)))


!**********************************************************************************************************************************

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QcQb)

!        print *, "Tree_Mqcqb=", colorfactor23*tree23QcQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- MqcH+MqcV ---

!=======================================================qc*s (Higgs)+ qc*s (Vector)================================================

!Higgs:- from i-combination  * -from Vector-contraction 

        do j=sumstartH, sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_Bc = dconjg(LLR(xtype,ytype,j))
          propac = dconjg(DenSmH23(j))

!*************************************************************t*s******************************************************************

!(p1-radiation)*(p1(Higgs)-radiation)

	tree23QcS =  tree23QcS + (4*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - muu5**2 - x2))*
     &                           (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                           Coupl23_At*Coupl23_Bt)/
     &                           (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                           (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p1-radiation)*(p2(Higgs)-radiation)

	tree23QcS =  tree23QcS + (-8*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &                           (-colorfactor23*((Coupl23_Ac*Coupl23_Bc*V2radc*
     &                           Coupl23_At*Coupl23_Bt)/
     &                           (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                           (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p2-radiation)*(p1(Higgs)-radiation)

	tree23QcS =  tree23QcS + (4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*
     &                           (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                           Coupl23_Dt*Coupl23_Ct)/
     &                           (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                           (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2-radiation)*(p2(Higgs)-radiation)

	tree23QcS =  tree23QcS - ((4*p1k1 + sqrtS**2*(1 - 2*muu1**2 - 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1))*
     &                           (-colorfactor23*((Coupl23_Ac*Coupl23_Bc*V2radc*
     &                           Coupl23_Dt*Coupl23_Ct)/
     &                           (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                           (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!*************************************************************u*s******************************************************************
!(p1-radiation)*(p1(Higgs)-radiation)

	tree23QcS =  tree23QcS + (4*p1k1 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - muu5**2 - x1))*
     &                           (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                           Coupl23_Au*Coupl23_Bu)/
     &                           (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))) 

!(p1-radiation)*(p2(Higgs)-radiation)

	tree23QcS =  tree23QcS + (-4*p1k1 - 8*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &                           (-colorfactor23*(Coupl23_Ac*Coupl23_Bc*V2radc*
     &                           Coupl23_Au*Coupl23_Bu)/
     &                           (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))

!(p2-radiation)*(p1(Higgs)-radiation)

	tree23QcS =  tree23QcS + (8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*
     &                           (colorfactor23*(Coupl23_Ac*Coupl23_Bc*V1radc*
     &                           Coupl23_Du*Coupl23_Cu)/
     &                           (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2-radiation)*(p2(Higgs)-radiation)

	tree23QcS =  tree23QcS - ((4*p1k2 + sqrtS**2*(1 - 2*muu1**2 - 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - x2))*
     &                           (-colorfactor23*((Coupl23_Ac*Coupl23_Bc*V2radc*
     &                           Coupl23_Du*Coupl23_Cu)/
     &                           (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!**********************************************************Dipoles( qc*s(Higgs) )**************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QcS = tree23QcS - 

!p2-splitting

!        tree23QcS = tree23QcS -  

!      endif
  
!**********************************************************************************************************************************

        enddo

        do j= sumstartV, sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))


!Vector:- from i-combination  * -from Vector-contraction * - from Propagator

!*************************************************************t*s******************************************************************
!(p1-radiation)*(p1(Vector)-radiation)

	tree23QcS =  tree23QcS + ((32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(20*p1k2 + sqrtS**2*(5 + muu3**2 + 3*muu4**2 - muu5**2 - 7*x1 - 4*x2)) - 
     &                            4*p1k2*sqrtS**2*(-3 + 6*muu1**2 - 6*muu2**2 - 3*muu3**2 - muu4**2 + muu5**2 + 10*x1 + 4*x2) + 
     &                            sqrtS**4*(-17 - 2*muu3**2 - muu3**4 - 6*muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + 18*muu5**2 + 2*muu3**2*muu5**2 + 2*muu4**2*muu5**2 - muu5**4 + 
     &                            4*muu1**2*(4 + 2*muu3**2 - muu4**2 - 4*muu5**2 - 5*x1) + 14*x1 + 2*muu3**2*x1 - 4*muu4**2*x1 - 2*muu5**2*x1 + x1**2 + 4*muu2**2*(6 + 3*muu4**2 - 6*muu5**2 - 5*x1 - 4*x2) + 
     &                            7*x2 - muu3**2*x2 - muu4**2*x2 + muu5**2*x2 + 5*x1*x2 + 2*x2**2))*
     &                            (colorfactor23*(Coupl23_Ac*V1radc*Coupl23_Bc*
     &                            Coupl23_At*Coupl23_Bt)/
     &                            (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))/2.

!(p1-radiation)*(p2(Vector)-radiation)

	tree23QcS =  tree23QcS - ((32*p1k2**2 + 4*p1k2*sqrtS**2*(5 + 2*muu1**2 - 2*muu2**2 + muu3**2 + 3*muu4**2 - 3*muu5**2 - 2*x1 - 4*x2) + 
     &                            sqrtS**4*(5 - 2*muu3**2 + muu3**4 + 10*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 6*x1 - 2*muu3**2*x1 - 4*muu4**2*x1 + 
     &                            2*muu5**2*x1 + x1**2 + 4*muu1**2*(3 + muu3**2 - 4*muu4**2 - 3*muu5**2 - 7*x1 - 5*x2) - 7*x2 - 3*muu3**2*x2 - 3*muu4**2*x2 + 3*muu5**2*x2 + 3*x1*x2 + 2*x2**2 + 
     &                            4*muu2**2*(-1 + muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2)) + 4*p1k1*(4*p1k2 + sqrtS**2*(3 + 4*muu1**2 - 4*muu2**2 - muu3**2 + 5*muu4**2 + muu5**2 + 3*x1 + 2*x2)))*
     &                            (-colorfactor23*((Coupl23_Ac*V2radc*Coupl23_Bc*
     &                            Coupl23_At*Coupl23_Bt)/
     &                            (propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))/2.

!(p2-radiation)*(p1(Vector)-radiation)

	tree23QcS =  tree23QcS - ((32*p1k1**2 + 4*p1k1*(4*p1k2 + sqrtS**2*(-5 + 2*muu1**2 - 2*muu2**2 - 3*muu3**2 - muu4**2 + 3*muu5**2 - 4*x1)) + 
     &                            sqrtS**2*(4*p1k2*(-3 + 4*muu1**2 - 4*muu2**2 - 5*muu3**2 + muu4**2 - muu5**2 - 4*x1 - 3*x2) + 
     &                            sqrtS**2*(5 + 10*muu3**2 + muu3**4 - 2*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 + 3*x1 + 3*muu3**2*x1 - 
     &                            muu4**2*x1 - 3*muu5**2*x1 + 2*x1**2 + 6*muu3**2*x2 - 4*muu4**2*x2 + 4*muu5**2*x2 + 7*x1*x2 + 7*x2**2 - 4*muu1**2*(1 + 2*muu3**2 - muu4**2 - muu5**2 + x2) - 
     &                            4*muu2**2*(-3 + 4*muu3**2 - muu4**2 + 3*muu5**2 + 4*x1 + 5*x2))))*
     &                            (colorfactor23*(Coupl23_Ac*V1radc*Coupl23_Bc*
     &                            Coupl23_Dt*Coupl23_Ct)/
     &                            (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))/2.

!(p2-radiation)*(p2(Vector)-radiation)

	tree23QcS =  tree23QcS + ((32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(20*p1k2 - sqrtS**2*(3 + 6*muu1**2 - 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + 4*x1)) - 
     &                            4*p1k2*sqrtS**2*(5 + 3*muu3**2 + muu4**2 - muu5**2 + 6*x1 + x2) + 
     &                            sqrtS**4*(-17 - 6*muu3**2 - muu3**4 - 2*muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + 18*muu5**2 + 2*muu3**2*muu5**2 + 2*muu4**2*muu5**2 - muu5**4 + 13*x1 + muu3**2*x1 + 5*muu4**2*x1 - 
     &                            muu5**2*x1 + 2*x1**2 + 4*muu1**2*(6 + 3*muu3**2 - 6*muu5**2 - x1 - 5*x2) + 24*x2 + 2*muu3**2*x2 + 4*muu4**2*x2 - 4*muu5**2*x2 - 3*x1*x2 - 5*x2**2 - 
     &                            4*muu2**2*(-4 + muu3**2 - 2*muu4**2 + 4*muu5**2 + 3*x1 + 5*x2)))*
     &                            (-colorfactor23*((Coupl23_Ac*V2radc*Coupl23_Bc*
     &                            Coupl23_Dt*Coupl23_Ct)/
     &                            (propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))/2.

!*************************************************************u*s******************************************************************
!(p1-radiation)*(p1(Vector)-radiation)

	tree23QcS =  tree23QcS - ((32*p1k1**2 + 32*p1k2**2 + 4*p1k1*(20*p1k2 + sqrtS**2*(3 - 6*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - 4*x1 - 10*x2)) + 
     &                            4*p1k2*sqrtS**2*(5 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 - 7*x2) - 
     &                            sqrtS**4*(17 + 6*muu3**2 + muu3**4 + 2*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 18*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 7*x1 + muu3**2*x1 + muu4**2*x1 - 
     &                            muu5**2*x1 - 2*x1**2 - 4*muu2**2*(6 + 3*muu3**2 - 6*muu5**2 - 4*x1 - 5*x2) - 14*x2 + 4*muu3**2*x2 - 2*muu4**2*x2 + 2*muu5**2*x2 - 5*x1*x2 - x2**2 + 
     &                            4*muu1**2*(-4 + muu3**2 - 2*muu4**2 + 4*muu5**2 + 5*x2)))*(colorfactor23*(Coupl23_Ac*
     &                            V1radc*Coupl23_Bc*Coupl23_Au*
     &                            Coupl23_Bu)/
     &                            (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))/2.

!(p1-radiation)*(p2(Vector)-radiation)

	tree23QcS =  tree23QcS + ((32*p1k1**2 + 4*p1k1*(4*p1k2 + sqrtS**2*(5 + 2*muu1**2 - 2*muu2**2 + 3*muu3**2 + muu4**2 - 3*muu5**2 - 4*x1 - 2*x2)) + 
     &                            sqrtS**2*(4*p1k2*(3 + 4*muu1**2 - 4*muu2**2 + 5*muu3**2 - muu4**2 + muu5**2 + 2*x1 + 3*x2) + 
     &                            sqrtS**2*(5 + 10*muu3**2 + muu3**4 - 2*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 7*x1 - 3*muu3**2*x1 - 
     &                            3*muu4**2*x1 + 3*muu5**2*x1 + 2*x1**2 - 6*x2 - 4*muu3**2*x2 - 2*muu4**2*x2 + 2*muu5**2*x2 + 3*x1*x2 + x2**2 + 4*muu2**2*(-1 - 2*muu3**2 + muu4**2 + muu5**2 + x1 + x2) - 
     &                            4*muu1**2*(-3 + 4*muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 7*x2))))*
     &                            (-colorfactor23*((Coupl23_Ac*V2radc*Coupl23_Bc*
     &                            Coupl23_Au*Coupl23_Bu)/
     &                            (propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))/2.

!(p2-radiation)*(p1(Vector)-radiation)

	tree23QcS =  tree23QcS + ((32*p1k2**2 + 4*p1k1*(4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 - 4*muu2**2 + muu3**2 - 5*muu4**2 - muu5**2 - 3*x1 - 4*x2)) + 
     &                            4*p1k2*sqrtS**2*(-5 + 2*muu1**2 - 2*muu2**2 - muu3**2 - 3*muu4**2 + 3*muu5**2 - 4*x2) + 
     &                            sqrtS**4*(5 - 2*muu3**2 + muu3**4 + 10*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 + 
     &                            4*muu1**2*(-1 + muu3**2 - 2*muu4**2 + muu5**2 - x1) - 4*muu3**2*x1 + 6*muu4**2*x1 + 4*muu5**2*x1 + 7*x1**2 + 4*muu2**2*(3 + muu3**2 - 4*muu4**2 - 3*muu5**2 - 5*x1 - 4*x2) + 
     &                            3*x2 - muu3**2*x2 + 3*muu4**2*x2 - 3*muu5**2*x2 + 7*x1*x2 + 2*x2**2))*
     &                            (colorfactor23*(Coupl23_Ac*V1radc*Coupl23_Bc*
     &                            Coupl23_Du*Coupl23_Cu)/
     &                            (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2))))/2.

!(p2-radiation)*(p2(Vector)-radiation)

	tree23QcS =  tree23QcS - ((32*p1k1**2 + 32*p1k2**2 - 4*p1k2*sqrtS**2*(3 + 6*muu1**2 - 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 + 4*x2) + 
     &                            sqrtS**4*(-17 - 2*muu3**2 - muu3**4 - 6*muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + 18*muu5**2 + 2*muu3**2*muu5**2 + 2*muu4**2*muu5**2 - muu5**4 + 24*x1 + 4*muu3**2*x1 + 
     &                            2*muu4**2*x1 - 4*muu5**2*x1 - 5*x1**2 + 4*muu2**2*(4 + 2*muu3**2 - muu4**2 - 4*muu5**2 - 5*x1 - 3*x2) + 4*muu1**2*(6 + 3*muu4**2 - 6*muu5**2 - 5*x1 - x2) + 13*x2 + 
     &                            5*muu3**2*x2 + muu4**2*x2 - muu5**2*x2 - 3*x1*x2 + 2*x2**2) + 4*p1k1*(20*p1k2 - sqrtS**2*(5 + muu3**2 + 3*muu4**2 - muu5**2 + x1 + 6*x2)))*
     &                            (-colorfactor23*((Coupl23_Ac*V2radc*Coupl23_Bc*
     &                            Coupl23_Du*Coupl23_Cu)/
     &                            (propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))/2.

!**********************************************************Dipoles( qc*s(Vector) )*************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QcS = tree23QcS - 

!p2-splitting

!        tree23QcS = tree23QcS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
         
        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + dreal(tree23QcS)

!   	print *, "Tree_Mqcs=", tree23QcS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mqct+Mqcu ---

!====================================================qc*t + qc*u (sFermion)========================================================

        do j= 1,2       

          Coupl23_Atc = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl23_Dtc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmtc = MSf(j,proptypetc,propgentc)
          sfwtc = WSf(j,proptypetc,propgentc)

          Coupl23_Auc = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl23_Duc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmuc = MSf(j,proptypeuc,propgenuc)
          sfwuc = WSf(j,proptypeuc,propgenuc)

          V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
          V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))  
   

! - from i-combination * -from Vector-contraction 

!(p1(t)-radiation)*(p1(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 - x2))*
     &  ((V1radc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_At*
     &      Coupl23_Bt)/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p1(t)-radiation)*(p2(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(12*p1k1 + 8*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &  (-((V2radc*Coupl23_Atc*
     &        Coupl23_Dtc*Coupl23_At*
     &        Coupl23_Bt)/
     &      ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p1(t)-radiation)*(channel(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k1 - 8*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + x1 - 4*x2))*
     &  ((V3tradc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_At*
     &      Coupl23_Bt)/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!(p1(u)-radiation)*(p1(t)-radiation)

        tree23QcTU =  tree23QcTU + sqrtS**2*(8*p1k1 + 12*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x2))*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &  ((V1radc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Au*
     &      Coupl23_Bu)/
     &    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p1(u)-radiation)*(p2(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu3**2 + muu5**2 + x1 + x2))*(8*p1k1 + 4*p1k2 - sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &  (-((V2radc*Coupl23_Atc*
     &        Coupl23_Dtc*Coupl23_Au*
     &        Coupl23_Bu)/
     &      ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p1(u)-radiation)*(channel(t)-radiation)

        tree23QcTU =  tree23QcTU + sqrtS**2*(-8*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - x2))*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &  ((V3tradc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Au*
     &      Coupl23_Bu)/
     &    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!(p2(t)-radiation)*(p1(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*(8*p1k1 + 12*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 - x2))*
     &  ((V1radc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Dt*
     &      Coupl23_Ct)/
     &    ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(t)-radiation)*(p2(t)-radiation)

        tree23QcTU =  tree23QcTU + (-4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 + x1))*(-8*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*
     &  (-((V2radc*Coupl23_Atc*
     &        Coupl23_Dtc*Coupl23_Dt*
     &        Coupl23_Ct)/
     &      ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p2(t)-radiation)*(channel(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*(-8*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - x2))*
     &  ((V3tradc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Dt*
     &      Coupl23_Ct)/
     &    ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))

!(p2(u)-radiation)*(p1(t)-radiation)

        tree23QcTU =  tree23QcTU + (4*p1k1 + 8*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 - muu3**2 + muu4**2 + muu5**2 - x1 - 2*x2))*(4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu4**2 + muu5**2 - x1 - x2))*
     &      ((V1radc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Du*
     &      Coupl23_Cu)/
     &    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p2(u)-radiation)*(p2(t)-radiation)

        tree23QcTU =  tree23QcTU + sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(-12*p1k1 - 8*p1k2 + sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 + 5*x1 + 4*x2))*
     &  (-((V2radc*Coupl23_Atc*
     &        Coupl23_Dtc*Coupl23_Du*
     &        Coupl23_Cu)/
     &      ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p2(u)-radiation)*(channel(t)-radiation)

        tree23QcTU =  tree23QcTU + sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(-4*p1k1 + 8*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 + x1 - 4*x2))*
     &  ((V3tradc*Coupl23_Atc*
     &      Coupl23_Dtc*Coupl23_Du*
     &      Coupl23_Cu)/
     &    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &      (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))))


!**********************************************************Dipoles( qc*t )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QcTU = tree23QcTU - 

!p2-splitting

!        tree23QcTU = tree23QcTU -  

!      endif
  
!**********************************************************************************************************************************
 
! - from i-combination * -from Vector-contraction * - from Propagator 

!(p1(t)-radiation)*(p1(u)-radiation)
     
        tree23QcTU =  tree23QcTU + sqrtS**2*(12*p1k1 + 8*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - x1))*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*
     &  ((V1radc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_At*
     &      Coupl23_Bt)/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))

!(p1(t)-radiation)*(p2(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu4**2 + muu5**2 + x1 + x2))*(4*p1k1 + 8*p1k2 - sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &  (-((V2radc*Coupl23_Duc*
     &        Coupl23_Auc*Coupl23_At*
     &        Coupl23_Bt)/
     &      ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p1(t)-radiation)*(channel(u)-radiation)
     
        tree23QcTU =  tree23QcTU + sqrtS**2*(4*p1k1 - 8*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - x1))*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*
     &  ((V3uradc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_At*
     &      Coupl23_Bt)/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!(p2(t)-radiation)*(p1(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*(4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu3**2 + muu5**2 - x1 - x2))*
     &  ((V1radc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Dt*
     &      Coupl23_Ct)/
     &    ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))
 
!(p2(t)-radiation)*(p2(u)-radiation)
     
        tree23QcTU =  tree23QcTU  - (sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*(8*p1k1 + 12*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 + 4*x1 + 5*x2))*
     &    (-((V2radc*Coupl23_Duc*
     &          Coupl23_Auc*Coupl23_Dt*
     &          Coupl23_Ct)/
     &        ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &          (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2)))))

!(p2(t)-radiation)*(channel(u)-radiation)
     
        tree23QcTU =  tree23QcTU + sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*(8*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 + x2))*
     &  ((V3uradc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Dt*
     &      Coupl23_Ct)/
     &    ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!(p1(u)-radiation)*(p1(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 - x1))*(8*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 2*x1 - x2))*
     &  ((V1radc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Au*
     &      Coupl23_Bu)/
     &    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))
 
!(p1(u)-radiation)*(p2(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(8*p1k1 + 12*p1k2 - sqrtS**2*(-5 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + 3*muu5**2 + 4*x1 + 5*x2))*
     &  (-((V2radc*Coupl23_Duc*
     &        Coupl23_Auc*Coupl23_Au*
     &        Coupl23_Bu)/
     &      ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p1(u)-radiation)*(channel(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(-8*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 6*muu2**2 + 3*muu3**2 + muu4**2 - muu5**2 - 4*x1 + x2))*
     &  ((V3uradc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Au*
     &      Coupl23_Bu)/
     &    ((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!(p2(u)-radiation)*(p1(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*(12*p1k1 + 8*p1k2 + sqrtS**2*(-5 + 2*muu1**2 + 6*muu2**2 + muu3**2 - muu4**2 + 3*muu5**2 - x1))*
     &  ((V1radc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Du*
     &      Coupl23_Cu)/
     &    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 + m2**2 + m5**2 - m1**2)))
 
!(p2(u)-radiation)*(p2(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 + x2))*(-4*p1k1 - 8*p1k2 + sqrtS**2*(-3 + 6*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + x1 + 2*x2))*
     &  (-((V2radc*Coupl23_Duc*
     &        Coupl23_Auc*Coupl23_Du*
     &        Coupl23_Cu)/
     &      ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 + m2**2))))

!(p2(u)-radiation)*(channel(u)-radiation)
     
        tree23QcTU =  tree23QcTU + (4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*(4*p1k1 - 8*p1k2 + sqrtS**2*(-1 + 6*muu1**2 + 2*muu2**2 + muu3**2 + 3*muu4**2 - muu5**2 - x1))*
     &  ((V3uradc*Coupl23_Duc*
     &      Coupl23_Auc*Coupl23_Du*
     &      Coupl23_Cu)/
     &    ((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))))

!**********************************************************Dipoles( qc*u )*********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QcTU = tree23QcTU - 

!p2-splitting

!        tree23QcTU = tree23QcTU -  

!      endif
  
!**********************************************************************************************************************************
    
        enddo

        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QcTU)

!    	print *, "Tree_Mqct+Tree_Mqcu=", colorfactor23*tree23QcTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

!==================================================================================================================================

      enddo

! --- Mss --

!*************************************************************Higgs-Higgs (s-Kanal)************************************************

!Factor 4 from polarisation-vector-contraction

      tree23SS = 0d0

      do i= sumstartH, sumendH       

        Coupl23_A = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
        Coupl23_B = LLR(xtype,ytype,i)
        propa = DenSmH23(i)
                       
         do j= sumstartH, sumendH
           
           Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)) * onshell_H(j)
           Coupl23_Bc = dconjg(LLR(xtype,ytype,j))
           propac = dconjg(DenSmH23(j))

!-from Vector-contraction

!(p1-radiation)*(p1-radiation)

	tree23SS = tree23SS - 4*((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                        V1rad)*Coupl23_A*Coupl23_B*Coupl23_Ac
     &                       *Coupl23_Bc/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!Ghosts

        tree23SS = tree23SS + 2d0*((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc
     &                       *V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *GGH(ytype,xtype,j)/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))

!(p1-radiation)*(p2-radiation)

	tree23SS = tree23SS + 8*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                        V1rad)*Coupl23_A*Coupl23_B*Coupl23_Ac
     &                       *Coupl23_Bc/
     &                        (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 - 
     &                        m1**2 + m1**2)*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 - 
     &                        m2**2 + m2**2))) )

!Ghosts

!(p1-radiation)*(p2-radiation)

        tree23SS = tree23SS - 2d0*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc
     &                       *V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *GGH(ytype,xtype,j)/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p2-radiation)*(p1-radiation)

        tree23SS = tree23SS - 2d0*(((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2rad
     &                       *V1radc)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *GGH(ytype,xtype,j)/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p2-radiation)*(p2-radiation)

	tree23SS = tree23SS + 4*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_A*Coupl23_B*Coupl23_Ac
     &                       *Coupl23_Bc/
     &                        (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                        sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))

!Ghosts

!(p2-radiation)*(p2-radiation)

        tree23SS = tree23SS - 2d0*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc
     &                       *V2rad)
     &                       *Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *GGH(ytype,xtype,j)/(propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))



!**********************************************************Dipoles( Higgs*Higgs )**************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting(Vector-Ghost)

!        tree23SS = tree23SS - 

!p2-splitting(Vector-Ghost)

!        tree23SS = tree23SS -  

!      endif
  
!**********************************************************************************************************************************

        enddo  

!*********************************************************Higgs-Vector (s-Kanal)***************************************************

!Gluon-final-states zero wire 0-couplings, W+/W--exchange =0 wire couplings

        do j=sumstartV, sumendV

           Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
           Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
           propac = dconjg(DenSmV23(j))

!- from Vector-Propagator * -from Vector-contraction   

!(p1-radiation)*(p1-radiation)
 
        tree23SS = tree23SS + 2*((-3*(4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*
     &                        V1radc*V1rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac
     &                       *Coupl23_Bc/
     &                        (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))    )

!Ghosts

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

!(p1(Higgs)-radiation)*(p1(Vector)-radiation)

        tree23SS = tree23SS - 2d0*(  -((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*V1radc*
     &                        V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/
     &                        (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   ) 

!(p1(Vector)-radiation)*(p1(Higgs)-radiation)

        tree23SS = tree23SS - 2d0*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*V1rad*
     &                        V1radc)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/
     &                        (2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))    )

        else

!(p1(Higgs)-radiation)*(p1(Vector)-radiation)


        tree23SS = tree23SS + 2d0*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*V1radc*
     &                        V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/
     &                        (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  ) 

!(p1(Vector)-radiation)*(p1(Higgs)-radiation)

        tree23SS = tree23SS - 2d0*(  ( (4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*V1rad*
     &                        V1radc)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/
     &                        (2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))   )

        endif

!(p1(Higgs)-radiation)*(p2(Vector)-radiation)
 
        tree23SS = tree23SS + 2*((3*sqrtS**2*(2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        V2radc*V1rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac
     &                       *Coupl23_Bc/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!Ghosts

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

!(p1(Higgs)-radiation)*(p2(Vector)-radiation)


        tree23SS = tree23SS + 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*V2radc*
     &                        V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/
     &                        (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )   

!(p2(Vector)-radiation)*(p1(Higgs)-radiation)

        tree23SS = tree23SS - 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*V1radc*
     &                        V2rad)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   ) 

        else

!(p1(Higgs)-radiation)*(p2(Vector)-radiation)


        tree23SS = tree23SS + 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*V2radc*
     &                        V1rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/
     &                        (2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )   

!(p2(Vector)-radiation)*(p1(Higgs)-radiation)

        tree23SS = tree23SS - 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*V1radc*
     &                        V2rad)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))   ) 

        endif

!(p2(Higgs)-radiation)*(p1(Vector)-radiation)
 
        tree23SS = tree23SS + 2*((-3*sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(2*p1k1 - 2*p1k2 + sqrtS**2*(muu3**2 - muu4**2 - x1 + x2))*
     &                        V1radc*V2rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac
     &                       *Coupl23_Bc/(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2))     )

!Ghosts

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

!(p2(Higgs)-radiation)*(p1(Vector)-radiation)

        tree23SS = tree23SS - 2d0*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(-4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 + x1 - x2))*V1radc*
     &                        V2rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/(2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 
     &                        )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  ) 

!(p1(Vector)-radiation)*(p2(Higgs)-radiation)


        tree23SS = tree23SS + 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*V2radc*
     &                        V1rad)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )  

        else

!(p2(Higgs)-radiation)*(p1(Vector)-radiation)

        tree23SS = tree23SS - 2d0*(  (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(-4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 + x2))*V1radc*
     &                        V2rad)*Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/(2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 
     &                        )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  ) 

!(p1(Vector)-radiation)*(p2(Higgs)-radiation)


        tree23SS = tree23SS + 2d0*(  (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*V2radc*
     &                        V1rad)*dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )  

        endif

!(p2(Higgs)-radiation)*(p2(Vector)-radiation)
 
        tree23SS = tree23SS + 2*((3*(2*p1k1 - 2*p1k2 + (-muu3**2 + muu4**2)*sqrtS**2)*
     &                        (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_A*Coupl23_B*Coupl23_Ac
     &                       *Coupl23_Bc/
     &                        (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                        sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!Ghosts

        if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

!(p2(Higgs)-radiation)*(p2(Vector)-radiation)

        tree23SS = tree23SS - 2d0*( (  (4*p1k1 + sqrtS**2*(1 - muu3**2 + muu4**2 - muu5**2 - x1 - x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)
     &                       *Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/(2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )   

!(p2(Vector)-radiation)*(p2(Higgs)-radiation)


        tree23SS = tree23SS + 2d0*(  ((4*p1k2 - sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)
     &                       *dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   ) 

        else

!(p2(Higgs)-radiation)*(p2(Vector)-radiation)

        tree23SS = tree23SS - 2d0*( (  (4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)
     &                       *Coupl23_A*GGH(xtype,ytype,i)*Coupl23_Ac
     &                       *dconjg(GGV(xtype,ytype,j))/(2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )   

!(p2(Vector)-radiation)*(p2(Higgs)-radiation)


        tree23SS = tree23SS + 2d0*(  ((4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        V2radc*V2rad)
     &                       *dconjg(Coupl23_A)*GGH(ytype,xtype,i)*VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
     &                       *GGV(xtype,ytype,j)/(2.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   ) 

        endif

!**********************************************************Dipoles( Higgs*Vector )*************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting(Vector-Ghost)

!        tree23SS = tree23SS - 

!p2-splitting(Vector-Ghost)

!        tree23SS = tree23SS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
      enddo

!****************************************************Vector-Vector (s-Kanal)*******************************************************
    
      do i=sumstartV, sumendV

        Coupl23_A = VSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
        Coupl23_B = VVV(i,xtype,ytype)
        propa = DenSmV23(i)
      
         do j=sumstartV, sumendV
           
           Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j)
           Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
           propac = dconjg(DenSmV23(j))

!-from Vector-contraction
 
!(p1-radiation)*(p1-radiation)
!-from amplitude

        tree23SS = tree23SS + ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                        (16*p1k1**2 + 16*p1k2**2 + 4*p1k1*(28*p1k2 + sqrtS**2*(-1 - 7*muu3**2 + 3*muu4**2 + muu5**2 - 3*x1 - 13*x2)) + 
     &                        4*p1k2*sqrtS**2*(-1 + 3*muu3**2 - 7*muu4**2 + muu5**2 - 13*x1 - 3*x2) + 
     &                        sqrtS**4*(1 - 5*muu3**4 + 2*muu4**2 - 5*muu4**4 - 2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 8*muu4**2*x1 + 3*x1**2 + 
     &                        2*muu3**2*(1 + 5*muu4**2 - muu5**2 + 6*x1 - 4*x2) + 12*muu4**2*x2 + 26*x1*x2 + 3*x2**2 - 
     &                        8*muu2**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2))))*V1radc*
     &                        V1rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!Ghosts
!(p1-radiation)*(p1-radiation)

        tree23SS = tree23SS - 2d0*((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*
     &                        (4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*V1radc*V1rad)
     &                       *Coupl23_A*GGV(xtype,ytype,i)
     &                       *Coupl23_Ac*dconjg(GGV(xtype,ytype,j))
     &                       /(4.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))

!(p1-radiation)*(p2-radiation)

        tree23SS = tree23SS + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        (-16*p1k1**2 - 16*p1k2**2 + 4*p1k2*sqrtS**2*(7*x1 + 2*x2) + p1k1*(-112*p1k2 + 4*sqrtS**2*(2*x1 + 7*x2)) + 
     &                        sqrtS**4*(17 - 2*muu3**2 - 5*muu3**4 - 2*muu4**2 + 10*muu3**2*muu4**2 - 5*muu4**4 - 18*muu5**2 - 2*muu3**2*muu5**2 - 
     &                        2*muu4**2*muu5**2 + muu5**4 - 17*x1 + 5*muu3**2*x1 - 5*muu4**2*x1 + muu5**2*x1 - 17*x2 - 5*muu3**2*x2 + 5*muu4**2*x2 + 
     &                        muu5**2*x2 + 4*muu1**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2)) + 
     &                        4*muu2**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2))))*V2radc*
     &                        V1rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(2.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   ) 

!Ghosts

!(p1-radiation)*(p2-radiation)

        tree23SS = tree23SS - 2d0*( (sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*(4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*
     &                        (-4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 + x2))*V2radc*V1rad)
     &                       *Coupl23_A*GGV(xtype,ytype,i)
     &                       *Coupl23_Ac*dconjg(GGV(xtype,ytype,j))
     &                       /(4.*propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2-radiation)*(p1-radiation)

        tree23SS = tree23SS - 2d0*( (sqrtS**2*(-2*muu1**2 - 2*muu2**2 + muu5**2 + x1 + x2)*(-4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 + x2))*
     &                        (-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*V1radc*V2rad)
     &                       *dconjg(Coupl23_A)*onshell_V(j)*dconjg(GGV(xtype,ytype,i))
     &                       *VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)*GGV(xtype,ytype,j)
     &                       /(4.*dconjg(propa)*dconjg(propac)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2  ))  )

!(p2-radiation)*(p2-radiation)

        tree23SS = tree23SS - ((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        (16*p1k1**2 + 16*p1k2**2 + 4*p1k1*(28*p1k2 + sqrtS**2*(1 + 7*muu3**2 - 3*muu4**2 - muu5**2 - x1 - x2)) - 
     &                        4*p1k2*sqrtS**2*(-1 + 3*muu3**2 - 7*muu4**2 + muu5**2 + x1 + x2) + 
     &                        sqrtS**4*(1 - 5*muu3**4 + 2*muu4**2 - 5*muu4**4 - 2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 2*x1 - 2*muu4**2*x1 + 2*muu5**2*x1 + 
     &                        x1**2 + 2*muu3**2*(1 + 5*muu4**2 - muu5**2 - x1 - x2) - 2*x2 - 2*muu4**2*x2 + 2*muu5**2*x2 + 2*x1*x2 + x2**2 - 
     &                        8*muu1**2*(muu3**2 + muu4**2 + 4*(-1 + muu5**2 + x1 + x2))))*V2radc*
     &                        V2rad)*Coupl23_A*Coupl23_B
     &                       *Coupl23_Ac*Coupl23_Bc
     &                       /(2.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                        sqrtS**2*x2  )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))

!Ghosts
!(p2-radiation)*(p2-radiation)

        tree23SS = tree23SS + 2d0*((4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*(4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                        (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)
     &                       *Coupl23_A*GGV(xtype,ytype,i)
     &                       *Coupl23_Ac*dconjg(GGV(xtype,ytype,j))
     &                       /(4.*propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))



!**********************************************************Dipoles( Vector*Vector )************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting(Vector-Ghost)

!        tree23SS = tree23SS - 

!p2-splitting(Vector-Ghost)

!        tree23SS = tree23SS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
      enddo

!**********************************************************************************************************************************

      QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23SS)

!        print *, "Tree_Mss=", colorfactor23*tree23SS
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mtt+Muu+Mtu+Mut --

!==========================================================sFermion (t- und u-Kanal)===============================================


      do i= 1, 2

        tree23TU = 0d0    
        tree23QTU = 0d0
        tree23STU = 0d0

        Coupl23_At = VSfSft(xtype,i,sfmeigen1,proptypet,sftype1,propgent,sfgen1)
        Coupl23_Dt = VSfSfu(ytype,sfmeigen2,i,sftype2,proptypet,sfgen2,propgent)
        sfmt = MSf(i,proptypet,propgent)
        sfwt = WSf(i,proptypet,propgent)

        Coupl23_Au = VSfSfu(ytype,i,sfmeigen1,proptypeu,sftype1,propgenu,sfgen1)
        Coupl23_Du = VSfSft(xtype,sfmeigen2,i,sftype2,proptypeu,sfgen2,propgenu)
        sfmu = MSf(i,proptypeu,propgenu)
        sfwu = WSf(i,proptypeu,propgenu)

        V3trad = VSfSf(ztype,i,i,proptypet,proptypet,propgent,propgent)
        V3urad = VSfSf(ztype,i,i,proptypeu,proptypeu,propgenu,propgenu) 

        do j=1, 2        

          Coupl23_Atc = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl23_Dtc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmtc = MSf(j,proptypetc,propgentc)
          sfwtc = WSf(j,proptypetc,propgentc)

          Coupl23_Auc = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl23_Duc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmuc = MSf(j,proptypeuc,propgenuc)
          sfwuc = WSf(j,proptypeuc,propgenuc)

          V3tradc = dconjg(VSfSf(ztype,j,j,proptypetc,proptypetc,propgentc,propgentc))
          V3uradc = dconjg(VSfSf(ztype,j,j,proptypeuc,proptypeuc,propgenuc,propgenuc))

 
!***********************************************************t*t********************************************************************
! -from Vector-contraction

!(p1-radiation)*(p1-radiation)

        tree23TU = tree23TU - colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                        (4*p1k1 + 8*p1k2 + sqrtS**2*(4*muu2**2 - muu3**2 + 2*(-1 + muu4**2 + muu5**2 - x2)))*
     &                        (4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*V1radc*
     &                        V1rad)*Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Atc*Coupl23_Dtc
     &                       /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2-radiation)*(p2-radiation)

        tree23TU = tree23TU + colorfactor23*(  (((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(8*p1k1 + 4*p1k2 + 
     &                        sqrtS**2*(2 - 4*muu1**2 - 2*muu3**2 + muu4**2 - 2*muu5**2 - 2*x1 - 2*x2))*
     &                        (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                        V2rad)*Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Atc*Coupl23_Dtc
     &                       /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))  )

!(channel-radiation)*(channel-radiation)

        tree23TU = tree23TU - colorfactor23*(  ((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*
     &                       (4*p1k1 - 4*p1k2 + sqrtS**2*(-2*muu1**2 - 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + 2*x2))*
     &                       V3tradc*V3trad)*Coupl23_At
     &                      *Coupl23_Dt
     &                      *Coupl23_Atc*Coupl23_Dtc
     &                      /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                       (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                       (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc)))  )

!(p1-radiation)*(p2-radiation)

        tree23TU = tree23TU - 2*(colorfactor23*(sqrtS**2*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 + x1))*
     &                        (4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 - x2))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        V2radc*V1rad)*Coupl23_At
     &                       *Coupl23_Dt
     &                       *Coupl23_Atc*Coupl23_Dtc
     &                       /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                        (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  ))   )

!(p1-radiation)*(channel-radiation)

        tree23TU = tree23TU - 2*(colorfactor23*((4*p1k2 + sqrtS**2*(4*muu2**2 + muu4**2 - 2*x2))*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                       (4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 - x2))*V3tradc*
     &                       V1rad)*Coupl23_At*Coupl23_Dt
     &                      *Coupl23_Atc*Coupl23_Dtc
     &                      /((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                       (m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                       (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!(p2-radiation)*(channel-radiation)

        tree23TU = tree23TU - 2*(colorfactor23*((4*p1k1 - (4*muu1**2 + muu3**2)*sqrtS**2)*(4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                        (4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 + x1))*V3tradc*
     &                        V2rad)*Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Atc*Coupl23_Dtc
     &                       /((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m3**2 - 2*p1k1 + m1**2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                        (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmtc**2+cI*sfmtc*sfwtc))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   )
 
!**********************************************************Dipoles( t*t )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23TU = tree23TU - 

!p2-splitting

!        tree23TU = tree23TU -  

!      endif
  
!*********************************************************t*u+u*t******************************************************************

! -from Vector-contraction

!TU+UT   

!(p1(T)-radiation)*(p1(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu3**2 + muu5**2 - x1 - x2))*
     &    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu4**2 + muu5**2 - x1 - x2))*V1radc*
     &    V1rad)*Coupl23_At*Coupl23_Dt*Coupl23_Auc*Coupl23_Duc/
     &  ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))          )
 

!(p2(T)-radiation)*(p2(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(     -(((4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu3**2 + muu5**2 + x1 + x2))*
     &      (4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 - muu4**2 + muu5**2 + x1 + x2))*
     &      (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*V2rad)*
     &       Coupl23_At*Coupl23_Dt
     &      *Coupl23_Auc*Coupl23_Duc/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)*
     &      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)))       )

!(channel(T)-radiation)*(channel(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(   (sqrtS**6*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &                        (-2 + 2*muu1**2 + 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2)*V3uradc*V3trad)*
     &             Coupl23_At*Coupl23_Dt
     &            *Coupl23_Auc*Coupl23_Duc/
     &  ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))         )

!(p1(T)-radiation)*(p2(U)-radiation)

       tree23TU = tree23TU - 2*colorfactor23*(   -((sqrtS**6*(2*muu1**2 + 2*muu2**2 + muu3**2 - 2*muu4**2 - 2*muu5**2 - x1)*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &      (-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*V2radc*V1rad)*
     &                        Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &    ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )
     &      *(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &        sqrtS**2*x1 - sqrtS**2*x2 )))         )

!(p2(T)-radiation)*(p1(U)-radiation)

       tree23TU = tree23TU - 2*colorfactor23*(    -((sqrtS**6*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(2*muu1**2 + 2*muu2**2 - 2*muu3**2 + muu4**2 - 2*muu5**2 - x2)*
     &      (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                        V2rad)*Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)
     &      *(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &        sqrtS**2*x1 - sqrtS**2*x2 )))        )


!(p1(T)-radiation)*(channel(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(    (sqrtS**2*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                        (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu4**2 + muu5**2 - x1 - x2))*(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &                        V3uradc*V1rad)*
     &                        Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &                       ((m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                       (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                       (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                       (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))        )

!(channel(T)-radiation)*(p1(U)-radiation)

       tree23TU = tree23TU - 2*colorfactor23*(         (sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &    (4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 - muu3**2 + muu5**2 - x1 - x2))*
     &                        V1radc*V3trad)*
     &                        Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &  ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &    (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &    (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )
 
!(p2(T)-radiation)*(channel(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(  -((sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 - muu3**2 + x1)
     &                         *(-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &      (-4*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 4*muu1**2 - muu3**2 + muu5**2 + x1 + x2))*V3uradc*V2rad)*
     &                        Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))      )
 
!(channel(T)-radiation)*(p2(U)-radiation)

        tree23TU = tree23TU - 2*colorfactor23*(  -((sqrtS**2*(-4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))
     &                       *(-2 + 2*muu1**2 + 2*muu2**2 - muu4**2 + x2)*
     &      (-4*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 4*muu1**2 - muu4**2 + muu5**2 + x1 + x2))*V2radc*V3trad)*
     &                        Coupl23_At*Coupl23_Dt
     &                       *Coupl23_Auc*Coupl23_Duc/
     &    ((m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &      (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &      (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2)))          )

!**********************************************************Dipoles( t*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23TU = tree23TU - 

!p2-splitting

!        tree23TU = tree23TU -  

!      endif
  
!***********************************************************u*u********************************************************************

! -from Vector-contraction

!(p1-radiation)*(p1-radiation)

        tree23TU = tree23TU - colorfactor23*(  ((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*(4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*
     &                        (8*p1k1 + 4*p1k2 + sqrtS**2*(-2 + 4*muu2**2 + 2*muu3**2 - muu4**2 + 2*muu5**2 - 2*x1))*
     &                        V1radc*V1rad)*Coupl23_Au
     &                       *Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))   )

!(p2-radiation)*(p2-radiation)

        tree23TU = tree23TU + colorfactor23*(  (((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                        (4*p1k1 + 8*p1k2 + sqrtS**2*(-4*muu1**2 + muu3**2 - 2*(-1 + muu4**2 + muu5**2 + x1 + x2)))*
     &                        V2radc*V2rad)*Coupl23_Au
     &                       *Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))   )

!(channel-radiation)*(channel-radiation)

        tree23TU = tree23TU - colorfactor23*(  ((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*
     &                        (-4*p1k1 + 4*p1k2 + sqrtS**2*(-2*muu1**2 - 2*muu2**2 - 2*muu3**2 - 2*muu4**2 + muu5**2 + 2*x1))*
     &                        V3uradc*V3urad)*Coupl23_Au
     &                       *Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc)))  )

!(p1-radiation)*(p2-radiation)

        tree23TU = tree23TU + 2*(colorfactor23*((sqrtS**2*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 - x1))*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                        (-4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 + x2))*V2radc*
     &                        V1rad)*Coupl23_Au*Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                        sqrtS**2*x1 - sqrtS**2*x2  )))  )

!(p1-radiation)*(channel-radiation)

        tree23TU = tree23TU - 2*(colorfactor23*((4*p1k1 + sqrtS**2*(4*muu2**2 + muu3**2 - 2*x1))*(4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                        (4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu5**2 - x1))*V3uradc*
     &                        V1rad)*Coupl23_Au*Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                        (m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2-radiation)*(channel-radiation)

        tree23TU = tree23TU - 2*(colorfactor23*((4*p1k2 - (4*muu1**2 + muu4**2)*sqrtS**2)*(4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                        (4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu4**2 - muu5**2 + x2))*V3uradc*
     &                        V2rad)*Coupl23_Au*Coupl23_Du
     &                       *Coupl23_Auc*Coupl23_Duc
     &                       /((m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*(m4**2 - 2*p1k2 + m1**2 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmuc**2+cI*sfmuc*sfwuc))*
     &                        (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!**********************************************************Dipoles( u*u )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23TU = tree23TU - 

!p2-splitting

!        tree23TU = tree23TU -  

!      endif
  
!**********************************************************************************************************************************

        enddo
 
        QQ2xx_M2Tree23VV = QQ2xx_M2Tree23VV + dreal(tree23TU)

!    	print *, "Tree_Mtt+Tree_Muu+Tree_Mut=", tree23TU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV 

! --- Mqt+Mtq+Mqu+Muq ---

!====================================================t*q + u*q (sFermion)==========================================================   

        Coupl23_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2))    

! - from i-combination * -from Vector-contraction

!(p1(t)-radiation)*(p1(q)-radiation)

        tree23QTU =  tree23QTU + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                            (4*p1k1 + 12*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + 3*muu4**2 + muu5**2 - x1 - 5*x2))*
     &                            V1radc*V1rad)
     &                           *Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt/
     &                            (2.*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))    )

!(p2(t)-radiation)*(p2(q)-radiation)

        tree23QTU =  tree23QTU + 2*(((12*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 8*muu1**2 + 3*muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                            (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))
     &                           *V2radc*
     &                            V2rad)*Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt/
     &                            (2.*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))    )

!(channel(t)-radiation)*(p1(q)-radiation)

        tree23QTU =  tree23QTU + 2*(((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                            (-4*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                            V1radc*V3trad)
     &                           *Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt/
     &                            (2.*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2))    )

!(channel(t)-radiation)*(p2(q)-radiation)


        tree23QTU =  tree23QTU - 2*(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                            (4*p1k1 - 4*p1k2 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                            V2radc*V3trad)
     &                           *Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt/
     &                            (2.*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))    )

!(p1(t)-radiation)*(p2(q)-radiation)

        tree23QTU =  tree23QTU - 2*((sqrtS**2*(4*p1k1 + 12*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + 3*muu4**2 + muu5**2 - x1 - 5*x2))*
     &                           (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc
     &                           *V1rad
     &                            )*Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt
     &                           /(2.*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))    )

!(p2(t)-radiation)*(p1(q)-radiation)

        tree23QTU =  tree23QTU - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-12*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 8*muu1**2 + 3*muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_At*Coupl23_Dt/
     &                            (2.*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))    )



!**********************************************************Dipoles( t*q )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QTU = tree23QTU - 

!p2-splitting

!        tree23QTU = tree23QTU -  

!      endif
  
!**********************************************************************************************************************************
 
! - from i-combination * -from Vector-contraction 

!(p1(u)-radiation)*(p1(q)-radiation)
     
        tree23QTU =  tree23QTU + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                           (12*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + 3*muu3**2 - muu4**2 + muu5**2 - 5*x1 - x2))*
     &                           V1radc*V1rad)*Coupl23_Ac
     &                          *Coupl23_Au*Coupl23_Du/
     &                           (2.*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                           (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(p2(u)-radiation)*(p2(q)-radiation)
     
        tree23QTU =  tree23QTU + 2*(((4*p1k1 + 12*p1k2 - sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + 3*muu4**2 + muu5**2 + x1 + x2))*
     &                            (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                            V2rad)*Coupl23_Ac
     &                           *Coupl23_Au*Coupl23_Du/
     &                            (2.*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!(channel(u)-radiation)*(p1(q)-radiation)
     
        tree23QTU =  tree23QTU + 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                            (4*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - x1 + x2))*
     &                            V1radc*V3urad)*Coupl23_Ac
     &                           *Coupl23_Au*Coupl23_Du/
     &                            (2.*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))     )

!(channel(u)-radiation)*(p2(q)-radiation)
     
        tree23QTU =  tree23QTU - 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                            (-4*p1k1 + 4*p1k2 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - x1 + x2))*
     &                            V2radc*V3urad)*Coupl23_Ac
     &                           *Coupl23_Au*Coupl23_Du/
     &                            (2.*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))     )

!(p1(u)-radiation)*(p2(q)-radiation)
     
        tree23QTU =  tree23QTU - 2*((sqrtS**2*(12*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + 3*muu3**2 - muu4**2 + muu5**2 - 5*x1 - x2))*
     &                            (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad
     &                            )*Coupl23_Ac
     &                           *Coupl23_Au*Coupl23_Du
     &                           /(2.*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2) *(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))     )

!(p2(u)-radiation)*(p1(q)-radiation)
     
        tree23QTU =  tree23QTU - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-4*p1k1 - 12*p1k2 + sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + 3*muu4**2 + muu5**2 + x1 + x2))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_Au*Coupl23_Du/
     &                            (2.*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))     )


!**********************************************************Dipoles( u*q )**********************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QTU = tree23QTU - 

!p2-splitting

!        tree23QTU = tree23QTU -  

!      endif
  
!**********************************************************************************************************************************     

        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23QTU)

!    	print *, "Tree_Mqt+Tree_Mqu=", colorfactor23*tree23QTU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

! --- Mst+Mts+Msu+Mus ---

!====================================================t*s (Higgs)+ t*s (Vector) + u*s (Higgs)+ u*s (Vector)=========================

        do j=sumstartH, sumendH

           Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
           Coupl23_Bc = dconjg(LLR(xtype,ytype,j))
           propac = dconjg(DenSmH23(j))


!**************************************************************t*s(Higgs)**********************************************************

! =0 for Gluon-final-state wire couplings

! -from Vector-contraction

!(p1(t)-radiation)*(p1(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                            (4*p1k1 + 12*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + 3*muu4**2 + muu5**2 - x1 - 5*x2))*
     &                            V1radc*V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2(t)-radiation)*(p2(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((12*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 8*muu1**2 + 3*muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                            (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                            V2rad)*Coupl23_Ac*Coupl23_Bc
     &                           *Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(channel(t)-radiation)*(p1(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                            (-4*p1k1 + 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                            V1radc*V3trad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(channel(t)-radiation)*(p2(Higgs)-radiation)


	tree23STU =  tree23STU + 2*(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                            (4*p1k1 - 4*p1k2 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1 - x2))*
     &                            V2radc*V3trad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   )

!(p1(t)-radiation)*(p2(Higgs)-radiation)


	tree23STU =  tree23STU + 2*((sqrtS**2*(4*p1k1 + 12*p1k2 + sqrtS**2*(-1 + 8*muu2**2 - muu3**2 + 3*muu4**2 + muu5**2 - x1 - 5*x2))*
     &                            (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad
     &                            )*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/(2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p1(Higgs)-radiation)


	tree23STU =  tree23STU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-12*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 8*muu1**2 + 3*muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )


!**********************************************************Dipoles( t*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23STU = tree23STU - 

!p2-splitting

!        tree23STU = tree23STU -  

!      endif
  
!**************************************************************u*s(Higgs)**********************************************************
        
! =0 for Gluon-final-state wire couplings

!(p1(u)-radiation)*(p1(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                            (12*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + 3*muu3**2 - muu4**2 + muu5**2 - 5*x1 - x2))*
     &                            V1radc*V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2(u)-radiation)*(p2(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((4*p1k1 + 12*p1k2 - sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + 3*muu4**2 + muu5**2 + x1 + x2))*
     &                            (4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                            V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(channel(u)-radiation)*(p1(Higgs)-radiation)

	tree23STU =  tree23STU - 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                            (4*p1k1 - 4*p1k2 + sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - x1 + x2))*
     &                            V1radc*V3urad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(channel(u)-radiation)*(p2(Higgs)-radiation)


	tree23STU =  tree23STU + 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                            (-4*p1k1 + 4*p1k2 - sqrtS**2*(-3 + 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - x1 + x2))*
     &                            V2radc*V3urad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))   )

!(p1(u)-radiation)*(p2(Higgs)-radiation)


	tree23STU =  tree23STU + 2*((sqrtS**2*(12*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 8*muu2**2 + 3*muu3**2 - muu4**2 + muu5**2 - 5*x1 - x2))*
     &                            (2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*V1rad
     &                            )*Coupl23_Ac*Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/(2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(u)-radiation)*(p1(Higgs)-radiation)


	tree23STU =  tree23STU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-4*p1k1 - 12*p1k2 + sqrtS**2*(-1 + 8*muu1**2 - muu3**2 + 3*muu4**2 + muu5**2 + x1 + x2))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )

!**********************************************************Dipoles( U*Higgs )******************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23STU = tree23STU - 

!p2-splitting

!        tree23STU = tree23STU -  

!      endif
  
!**********************************************************************************************************************************

        enddo

        do j=sumstartV, sumendV        

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl23_Bc = dconjg(VVV(j,xtype,ytype))
          propac = dconjg(DenSmV23(j))      

!**************************************************************t*s(Vector)*********************************************************     

!Vector: -(Propagator) * -from Vector-contraction  

!(p1(t)-radiation)*(p1(Vector)-radiation)

	tree23STU =  tree23STU + 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                            (16*p1k1**2 + 16*p1k2**2 + 8*p1k1*(2*p1k2 - sqrtS**2*(1 + muu3**2 + muu4**2 - muu5**2 + x1)) - 
     &                            8*p1k2*sqrtS**2*(2*muu3**2 + x1 + 2*x2) + sqrtS**4*
     &                            (1 + muu3**4 + 2*muu4**2 + muu4**4 - 2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 + 2*x1 + 2*muu4**2*x1 - 2*muu5**2*x1 + x1**2 - 
     &                            2*muu3**2*(-1 + muu4**2 + muu5**2 - x1 - 3*x2) - 2*x2 - 2*muu4**2*x2 + 2*muu5**2*x2 + 2*x1*x2 + 5*x2**2 - 
     &                            8*muu2**2*(-1 + muu3**2 + muu4**2 + muu5**2 + x1 + x2)))*V1radc*
     &                            V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2(t)-radiation)*(p2(Vector)-radiation)

	tree23STU =  tree23STU + 2*( ((-4*p1k1 - 4*p1k2 + sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                            (16*p1k1**2 + 16*p1k2**2 + 16*p1k1*(p1k2 + m4**2) + 8*p1k2*sqrtS**2*(1 + muu3**2 + muu4**2 - muu5**2 - x1 - x2) + 
     &                            sqrtS**4*(muu3**4 + (-1 - muu4**2 + muu5**2 + x1 + x2)**2 - 2*muu3**2*(-1 + muu4**2 + muu5**2 + x1 + x2) - 
     &                            8*muu1**2*(-1 + muu3**2 + muu4**2 + muu5**2 + x1 + x2)))*V2radc*
     &                            V2rad)*Coupl23_Ac*Coupl23_Bc
     &                           *Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(channel(t)-radiation)*(p1(Vector)-radiation)


	tree23STU =  tree23STU - 2*(((4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 - x2))*
     &                            (16*p1k2**2 - 4*p1k2*sqrtS**2*(-5 + 4*muu1**2 - 4*muu2**2 - muu3**2 - muu4**2 + muu5**2 + 4*x1 + 3*x2) + 
     &                            sqrtS**4*(-4 - 8*muu3**2 + 4*muu4**2 + 4*muu5**2 + x1 + muu3**2*x1 - 3*muu4**2*x1 - muu5**2*x1 + x1**2 + 
     &                            4*muu2**2*(3 + 2*muu3**2 + muu4**2 - 3*muu5**2 - 2*x1 - x2) - 6*x2 - 2*muu4**2*x2 + 2*muu5**2*x2 + 5*x1*x2 + 2*x2**2 + 
     &                            4*muu1**2*(1 + 2*muu3**2 - muu4**2 - muu5**2 + x2)) + 
     &                            4*p1k1*(8*p1k2 - sqrtS**2*(-2 + 2*muu1**2 - 2*muu2**2 - 2*muu4**2 + x1 + 3*x2)))*V1radc*
     &                            V3trad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(channel(t)-radiation)*(p2(Vector)-radiation)


	tree23STU =  tree23STU - 2*(((4*p1k1 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 + x1))*
     &                            (16*p1k1**2 + 4*p1k1*(8*p1k2 - sqrtS**2*(5 + 4*muu1**2 - 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x1)) + 
     &                            sqrtS**2*(-4*p1k2*(2 + 2*muu1**2 - 2*muu2**2 + 2*muu3**2 + x1 - x2) + 
     &                            sqrtS**2*(-4 + 4*muu3**2 - 8*muu4**2 + 4*muu5**2 + 4*x1 + 2*muu4**2*x1 + 
     &                            4*muu1**2*(3 + muu3**2 + 2*muu4**2 - 3*muu5**2 + x1 - x2) + 5*x2 + muu3**2*x2 + muu4**2*x2 - muu5**2*x2 - x1*x2 - x2**2 - 
     &                            4*muu2**2*(-1 + muu3**2 - 2*muu4**2 + muu5**2 + x1 + x2))))*V2radc*
     &                            V3trad)*Coupl23_Ac*Coupl23_Bc
     &                           *Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )

!(p1(t)-radiation)*(p2(Vector)-radiation)

	tree23STU =  tree23STU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-32*p1k2**2 + 4*p1k2*sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 - muu3**2 + muu4**2 - muu5**2 + 3*x1 + 4*x2) + 
     &                            sqrtS**4*(5 + 2*muu3**2 + muu3**4 - 6*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 2*muu4**2*muu5**2 + 
     &                            muu5**4 - 7*x1 + muu3**2*x1 - muu4**2*x1 - muu5**2*x1 + 4*muu1**2*(-1 + 2*muu4**2 + muu5**2 + 2*x1) + 
     &                            4*muu2**2*(-1 - 2*muu3**2 + 4*muu4**2 + muu5**2 + 2*x1) - 3*x2 + 3*muu3**2*x2 - 3*muu4**2*x2 + 3*muu5**2*x2) - 
     &                            4*p1k1*(4*p1k2 + sqrtS**2*(-1 + 2*muu1**2 - 2*muu2**2 + muu3**2 - muu4**2 - muu5**2 + x2)))*
     &                            V2radc*V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m4**2 + 2*p1k2 + m2**2 - sqrtS**2*x2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )

!(p2(t)-radiation)*(p1(Vector)-radiation)

	tree23STU =  tree23STU + 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-32*p1k1**2 - 4*p1k1*(4*p1k2 - sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 + 4*x1 - x2)) + 
     &                            sqrtS**2*(-4*p1k2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 - 3*x1) + 
     &                            sqrtS**2*(5 - 6*muu3**2 + muu3**4 + 2*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 
     &                            2*muu4**2*muu5**2 + muu5**4 - 5*x1 - muu3**2*x1 + muu4**2*x1 + muu5**2*x1 - 5*x2 + muu3**2*x2 - muu4**2*x2 + muu5**2*x2 + 
     &                            4*muu2**2*(-1 + 2*muu3**2 + muu5**2 + x1 + x2) + 4*muu1**2*(-1 + 4*muu3**2 - 2*muu4**2 + muu5**2 - x1 + 3*x2))))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_At
     &                           *Coupl23_Dt/
     &                            (2.*propac*(m3**2 - 2*p1k1 + m1**2 - (sfmt**2-cI*sfmt*sfwt))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )

!**********************************************************Dipoles( t*Vector )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23STU = tree23STU - 

!p2-splitting

!        tree23STU = tree23STU -  

!      endif
  
!**************************************************************u*s(Vector)*********************************************************

!(p1(u)-radiation)*(p1(Vector)-radiation)
         
	tree23STU =  tree23STU - 2*(((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*
     &                            (16*p1k1**2 + 16*p1k2**2 - 8*p1k2*sqrtS**2*(1 + muu3**2 + muu4**2 - muu5**2 + x2) + 
     &                            sqrtS**4*(1 + muu3**4 + 2*muu4**2 + muu4**4 - 2*muu5**2 - 2*muu4**2*muu5**2 + muu5**4 - 2*x1 + 6*muu4**2*x1 + 2*muu5**2*x1 + 
     &                            5*x1**2 - 2*muu3**2*(-1 + muu4**2 + muu5**2 + x1 - x2) + 2*x2 + 2*muu4**2*x2 - 2*muu5**2*x2 + 2*x1*x2 + x2**2 - 
     &                            8*muu2**2*(-1 + muu3**2 + muu4**2 + muu5**2 + x1 + x2)) + 8*p1k1*(2*p1k2 - sqrtS**2*(2*muu4**2 + 2*x1 + x2)))*
     &                            V1radc*V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(p2(u)-radiation)*(p2(Vector)-radiation)
         
	tree23STU =  tree23STU + 2*(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*
     &                            (16*p1k1**2 + 16*p1k2**2 + 16*muu3**2*p1k2*sqrtS**2 + 8*p1k1*(2*p1k2 + sqrtS**2*(1 + muu3**2 + muu4**2 - muu5**2 - x1 - x2)) + 
     &                            sqrtS**4*(muu3**4 + (-1 - muu4**2 + muu5**2 + x1 + x2)**2 - 2*muu3**2*(-1 + muu4**2 + muu5**2 + x1 + x2) - 
     &                            8*muu1**2*(-1 + muu3**2 + muu4**2 + muu5**2 + x1 + x2)))*V2radc*
     &                            V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )


!(channel(u)-radiation)*(p1(Vector)-radiation)
         
	tree23STU =  tree23STU + 2*(((4*p1k1 + sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 + muu3**2 - muu4**2 - x1))*
     &                            (16*p1k1**2 + 4*p1k1*(8*p1k2 + sqrtS**2*(5 - 4*muu1**2 + 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 - 3*x1 - 4*x2)) + 
     &                            sqrtS**2*(-4*p1k2*(-2 + 2*muu1**2 - 2*muu2**2 - 2*muu3**2 + 3*x1 + x2) + 
     &                            sqrtS**2*(-4 + 4*muu3**2 - 8*muu4**2 + 4*muu5**2 - 4*muu1**2*(-1 + muu3**2 - 2*muu4**2 + muu5**2 - x1) - 6*x1 - 2*muu3**2*x1 + 
     &                            2*muu5**2*x1 + 2*x1**2 + 4*muu2**2*(3 + muu3**2 + 2*muu4**2 - 3*muu5**2 - x1 - 2*x2) + x2 - 3*muu3**2*x2 + muu4**2*x2 - 
     &                            muu5**2*x2 + 5*x1*x2 + x2**2)))*V1radc*V3urad)
     &                           *Coupl23_Ac*Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))  )

!(channel(u)-radiation)*(p2(Vector)-radiation)
         
	tree23STU =  tree23STU + 2*(((4*p1k2 - sqrtS**2*(-1 + 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + x2))*
     &                            (16*p1k2**2 - 4*p1k2*sqrtS**2*(5 + 4*muu1**2 - 4*muu2**2 + muu3**2 + muu4**2 - muu5**2 + x2) + 
     &                            4*p1k1*(8*p1k2 - sqrtS**2*(2 + 2*muu1**2 - 2*muu2**2 + 2*muu4**2 - x1 + x2)) + 
     &                            sqrtS**4*(-4 - 8*muu3**2 + 4*muu4**2 + 4*muu5**2 + 5*x1 + muu3**2*x1 + muu4**2*x1 - muu5**2*x1 - x1**2 + 4*x2 + 2*muu3**2*x2 - 
     &                            x1*x2 + 4*muu1**2*(3 + 2*muu3**2 + muu4**2 - 3*muu5**2 - x1 + x2) + muu2**2*(8*muu3**2 - 4*(-1 + muu4**2 + muu5**2 + x1 + x2)))
     &                            )*V2radc*V3urad)
     &                           *Coupl23_Ac*Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 ))  )


!(p1(u)-radiation)*(p2(Vector)-radiation)

         
	tree23STU =  tree23STU - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-32*p1k1**2 - 4*p1k1*(4*p1k2 + sqrtS**2*(1 - 2*muu1**2 + 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 - 4*x1 - 3*x2)) + 
     &                            sqrtS**2*(-4*p1k2*(-1 + 2*muu1**2 - 2*muu2**2 - muu3**2 + muu4**2 - muu5**2 + x1) + 
     &                            sqrtS**2*(5 - 6*muu3**2 + muu3**4 + 2*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 
     &                            2*muu4**2*muu5**2 + muu5**4 - 3*x1 - 3*muu3**2*x1 + 3*muu4**2*x1 + 3*muu5**2*x1 - 7*x2 - muu3**2*x2 + muu4**2*x2 - 
     &                            muu5**2*x2 + 4*muu1**2*(-1 + 2*muu3**2 + muu5**2 + 2*x2) + 4*muu2**2*(-1 + 4*muu3**2 - 2*muu4**2 + muu5**2 + 2*x2))))*
     &                            V2radc*V1rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m3**2 + 2*p1k1 + m2**2 - sqrtS**2*x1 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )


!(p2(u)-radiation)*(p1(Vector)-radiation)

         
	tree23STU =  tree23STU - 2*((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*
     &                            (-32*p1k2**2 - 4*p1k1*(4*p1k2 + sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 - muu3**2 + muu4**2 + muu5**2 - 3*x2)) + 
     &                            4*p1k2*sqrtS**2*(1 + 2*muu1**2 - 2*muu2**2 + muu3**2 - muu4**2 + muu5**2 - x1 + 4*x2) + 
     &                            sqrtS**4*(5 + 2*muu3**2 + muu3**4 - 6*muu4**2 - 2*muu3**2*muu4**2 + muu4**4 - 6*muu5**2 - 2*muu3**2*muu5**2 - 
     &                            2*muu4**2*muu5**2 + muu5**4 - 5*x1 - muu3**2*x1 + muu4**2*x1 + muu5**2*x1 - 5*x2 + muu3**2*x2 - muu4**2*x2 + muu5**2*x2 - 
     &                            4*muu1**2*(1 + 2*muu3**2 - 4*muu4**2 - muu5**2 - 3*x1 + x2) + 4*muu2**2*(-1 + 2*muu4**2 + muu5**2 + x1 + x2)))*
     &                            V1radc*V2rad)*Coupl23_Ac
     &                           *Coupl23_Bc*Coupl23_Au
     &                           *Coupl23_Du/
     &                            (2.*propac*(m4**2 - 2*p1k2 + m1**2 - (sfmu**2-cI*sfmu*sfwu))*
     &                            (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - 
     &                            sqrtS**2*x1 - sqrtS**2*x2  ))  )


!**********************************************************Dipoles( t*Vector )*****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23STU = tree23STU - 

!p2-splitting

!        tree23STU = tree23STU -  

!      endif
  
!**********************************************************************************************************************************
         
        enddo

        QQ2xx_M2Tree23VV= QQ2xx_M2Tree23VV + colorfactor23*dreal(tree23STU)

!    	print *, "Tree_Mst+Tree_Msu=", colorfactor23*tree23STU
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV

!================================================================================================================================== 

      enddo         
  
!==================================================================================================================================        
 
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23VV
!       QQ2xx_M2Tree23VV= tree23SS   !Nur bestimmte Kanaele sollen zaehlen

       end
