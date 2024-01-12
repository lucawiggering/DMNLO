!**********************************************************************************************************************************

	double precision function QQ2xx_M2Tree23FFbar()

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to3.h"



	integer i,j
	complex*16 Coupl23_A,Coupl23_Ac,Coupl23_BL,Coupl23_BR,Coupl23_BLc,Coupl23_BRc,propa,propac,treeSS23,tree23QbQb,tree23QbS

        QQ2xx_M2Tree23FFbar = 0d0

!==================================================================================================================================
!           ---               Fermion-Antifermion-final-states  (for el-we final states, Gluon/Photon-radiation)               ---
! =================================================================================================================================

! --- qbqb --

!==================================================================qb*qb===========================================================
  
      do i=sumstartV ,sumendV

        tree23QbQb = 0d0
        tree23QbS = 0d0

        Coupl23_A = VVSfSf(ztype,i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)
        Coupl23_BL = BL(i,f4type,f3type,f4gen,f3gen)
        Coupl23_BR = BR(i,f4type,f3type,f4gen,f3gen)
        propa = DenSmV23(i)

        do j=sumstartV ,sumendV

          Coupl23_Ac = dconjg(VVSfSf(ztype,j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)) * onshell_V(j) 
          Coupl23_BLc = dconjg(BL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(BR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmV23(j))
 
! -from Vector-contraction 

        tree23QbQb = tree23QbQb + (-4*(-(sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BL*Coupl23_BLc)/2. - 
     &                            2*m3*m4*Coupl23_BR*Coupl23_BLc - 2*m3*m4*Coupl23_BL*Coupl23_BRc - 
     &                            (sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BR*Coupl23_BRc)/2.)*
     &                            Coupl23_Ac*Coupl23_A)/(propac*propa)

        enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23FFbar = QQ2xx_M2Tree23FFbar + NC*colorfactor23*dreal(tree23QbQb)

!        print *, "Tree_Mqbqb=", NC*colorfactor23*tree23QbQb
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23FFbar

! --- MqbH+MqbV --


!=====================================================qb*s (Higgs)+ qb*s (Vector)==================================================


!************************************************************Higgs*****************************************************************

!Higgs:- from i-combination  * -from Vector-contraction *- from Vectorpropagator

        do j=sumstartH ,sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_BLc = dconjg(DL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(DR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmH23(j))

!(p1(Higgs)-radiation)

	tree23QbS =  tree23QbS - 2*((colorfactor23*sqrtS*(Coupl23_BR*(muu4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*Coupl23_BLc + 
     &                            muu3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*Coupl23_BRc) + 
     &                            Coupl23_BL*(muu3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1))*Coupl23_BLc + 
     &                            muu4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x2))*Coupl23_BRc))*Coupl23_Ac*
     &                            V1radc*Coupl23_A)/
     &                            (propac*propa*( 2*p1k1 + 2*p1k2 - sqrtS**2  + m2**2 + m5**2 - m1**2)))

!(p2(Higgs)-radiation)

	tree23QbS =  tree23QbS - 2*((colorfactor23*sqrtS*(Coupl23_BR*(muu4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*Coupl23_BLc + 
     &                            muu3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - 2*x2))*Coupl23_BRc) + 
     &                            Coupl23_BL*(muu3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - 2*x2))*Coupl23_BLc + 
     &                            muu4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + 2*x1 + x2))*Coupl23_BRc))*
     &                            Coupl23_Ac*V2radc*
     &                            Coupl23_A)/
     &                            (propac*propa*( - 2*p1k1 - 2*p1k2 - sqrtS**2 + m1**2 + m5**2 + sqrtS**2*x1 + sqrtS**2*x2 - 
     &                            m2**2)))

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
!Vector:-(i-combination) * -from Vector-contraction   

        do j=sumstartV ,sumendV

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j)  
          Coupl23_BLc = dconjg(BL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(BR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmV23(j))


!(p1(Vector)-radiation)

	tree23QbS =  tree23QbS + 2*((colorfactor23*Coupl23_Ac*V1radc*
     &                            (Coupl23_BR*(muu3*muu4*sqrtS**2*(-8*p1k1 - 8*p1k2 + sqrtS**2*(4 - 2*muu1**2 - 6*muu2**2 - 2*muu5**2 + x1 + x2))*Coupl23_BLc + 
     &                            Coupl23_BRc*(16*p1k1*p1k2 + 8*muu4**2*p1k1*sqrtS**2 + 8*muu3**2*p1k2*sqrtS**2 - sqrtS**4 + muu1**2*sqrtS**4 + 3*muu2**2*sqrtS**4 - 
     &                            2*muu3**2*sqrtS**4 + muu1**2*muu3**2*sqrtS**4 + 3*muu2**2*muu3**2*sqrtS**4 - muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + muu1**2*muu4**2*sqrtS**4 + 
     &                            3*muu2**2*muu4**2*sqrtS**4 + 2*muu3**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + muu5**2*sqrtS**4 - muu1**2*muu5**2*sqrtS**4 - 3*muu2**2*muu5**2*sqrtS**4 + 
     &                            muu3**2*muu5**2*sqrtS**4 + muu4**2*muu5**2*sqrtS**4 - 6*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - muu1**2*sqrtS**4*x1 - 3*muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 - 
     &                            2*muu4**2*sqrtS**4*x1 - 6*p1k1*sqrtS**2*x2 + sqrtS**4*x2 - muu1**2*sqrtS**4*x2 - 3*muu2**2*sqrtS**4*x2 - 2*muu3**2*sqrtS**4*x2 + muu4**2*sqrtS**4*x2 + 2*sqrtS**4*x1*x2 )) + 
     &                            Coupl23_BL*(muu3*muu4*sqrtS**2*(-8*p1k1 - 8*p1k2 + sqrtS**2*(4 - 2*muu1**2 - 6*muu2**2 - 2*muu5**2 + x1 + x2))*Coupl23_BRc + 
     &                            Coupl23_BLc*(16*p1k1*p1k2 + 8*muu4**2*p1k1*sqrtS**2 + 8*muu3**2*p1k2*sqrtS**2 - sqrtS**4 + muu1**2*sqrtS**4 + 3*muu2**2*sqrtS**4 - 
     &                            2*muu3**2*sqrtS**4 + muu1**2*muu3**2*sqrtS**4 + 3*muu2**2*muu3**2*sqrtS**4 - muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + muu1**2*muu4**2*sqrtS**4 + 
     &                            3*muu2**2*muu4**2*sqrtS**4 + 2*muu3**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + muu5**2*sqrtS**4 - muu1**2*muu5**2*sqrtS**4 - 3*muu2**2*muu5**2*sqrtS**4 + 
     &                            muu3**2*muu5**2*sqrtS**4 + muu4**2*muu5**2*sqrtS**4 - 6*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - muu1**2*sqrtS**4*x1 - 3*muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 - 
     &                            2*muu4**2*sqrtS**4*x1 - 6*p1k1*sqrtS**2*x2 + sqrtS**4*x2 - muu1**2*sqrtS**4*x2 - 3*muu2**2*sqrtS**4*x2 - 2*muu3**2*sqrtS**4*x2 + muu4**2*sqrtS**4*x2 + 2*sqrtS**4*x1*x2 )))
     &                            *Coupl23_A)/                                                                                                                                                  
     &                            (propac*propa*( 2*p1k1 + 2*p1k2 - sqrtS**2  + m2**2 + m5**2 - m1**2)))

!(p2(Vector)-radiation)

	tree23QbS =  tree23QbS +  2*((colorfactor23*Coupl23_Ac*V2radc*
     &                             (Coupl23_BL*(-(muu3*muu4*sqrtS**2*(-8*p1k1 - 8*p1k2 + sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                             Coupl23_BRc) + Coupl23_BLc*
     &                             (16*p1k1*p1k2 - 8*muu4**2*p1k1*sqrtS**2 - 8*muu3**2*p1k2*sqrtS**2 - sqrtS**4 + 3*muu1**2*sqrtS**4 + muu2**2*sqrtS**4 - 2*muu3**2*sqrtS**4 + 3*muu1**2*muu3**2*sqrtS**4 + 
     &                             muu2**2*muu3**2*sqrtS**4 - muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + 3*muu1**2*muu4**2*sqrtS**4 + muu2**2*muu4**2*sqrtS**4 + 2*muu3**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + 
     &                             muu5**2*sqrtS**4 - 3*muu1**2*muu5**2*sqrtS**4 - muu2**2*muu5**2*sqrtS**4 + muu3**2*muu5**2*sqrtS**4 + muu4**2*muu5**2*sqrtS**4 - 2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - 
     &                             3*muu1**2*sqrtS**4*x1 - muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + 2*muu4**2*sqrtS**4*x1 - 2*p1k1*sqrtS**2*x2 + sqrtS**4*x2 - 3*muu1**2*sqrtS**4*x2 - 
     &                             muu2**2*sqrtS**4*x2 + 2*muu3**2*sqrtS**4*x2 + muu4**2*sqrtS**4*x2 )) + 
     &                             Coupl23_BR*(-(muu3*muu4*sqrtS**2*(-8*p1k1 - 8*p1k2 + sqrtS**2*(-4 + 6*muu1**2 + 2*muu2**2 + 2*muu5**2 + 3*x1 + 3*x2))*
     &                             Coupl23_BLc) + Coupl23_BRc*
     &                             (16*p1k1*p1k2 - 8*muu4**2*p1k1*sqrtS**2 - 8*muu3**2*p1k2*sqrtS**2 - sqrtS**4 + 3*muu1**2*sqrtS**4 + muu2**2*sqrtS**4 - 2*muu3**2*sqrtS**4 + 3*muu1**2*muu3**2*sqrtS**4 + 
     &                             muu2**2*muu3**2*sqrtS**4 - muu3**4*sqrtS**4 - 2*muu4**2*sqrtS**4 + 3*muu1**2*muu4**2*sqrtS**4 + muu2**2*muu4**2*sqrtS**4 + 2*muu3**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 + 
     &                             muu5**2*sqrtS**4 - 3*muu1**2*muu5**2*sqrtS**4 - muu2**2*muu5**2*sqrtS**4 + muu3**2*muu5**2*sqrtS**4 + muu4**2*muu5**2*sqrtS**4 - 2*p1k2*sqrtS**2*x1 + sqrtS**4*x1 - 
     &                             3*muu1**2*sqrtS**4*x1 - muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 + 2*muu4**2*sqrtS**4*x1 - 2*p1k1*sqrtS**2*x2 + sqrtS**4*x2 - 3*muu1**2*sqrtS**4*x2 - 
     &                             muu2**2*sqrtS**4*x2 + 2*muu3**2*sqrtS**4*x2 + muu4**2*sqrtS**4*x2 )))*
     &                             Coupl23_A)/
     &                             (propac*propa*( - 2*p1k1 - 2*p1k2 - sqrtS**2 + m1**2 + m5**2 + sqrtS**2*x1 + sqrtS**2*x2 - 
     &                             m2**2)))

!**********************************************************Dipoles( Vector*qb )****************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        tree23QbS = tree23QbS - 

!p2-splitting

!        tree23QbS = tree23QbS -  

!      endif
  
!**********************************************************************************************************************************

        enddo
         
        QQ2xx_M2Tree23FFbar= QQ2xx_M2Tree23FFbar + NC*dreal(tree23QbS)

!    	print *, "Tree_Mqbs=", NC*tree23QbS
!     	print *, "Tree_M2Tree=", QQ2xx_M2Tree23FFbar


!==================================================================================================================================

      enddo

! --- Mss --

!==========================================================s-channel===============================================================

      treeSS23 = 0d0

!*******************************************************Higgs-exchange*************************************************************

!- from Polarisation-Vectors

      do i= sumstartH ,sumendH  

        Coupl23_A = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
        Coupl23_BL = DL(i,f4type,f3type,f4gen,f3gen)
        Coupl23_BR = DR(i,f4type,f3type,f4gen,f3gen)
        propa = DenSmH23(i)      
                          
        do j= sumstartH ,sumendH 

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)) * onshell_H(j)
          Coupl23_BLc = dconjg(DL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(DR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmH23(j))

!(p1-radiation)*(p1-radiation)

	treeSS23 = treeSS23 - (((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad))*(sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BLc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BRc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BLc*Coupl23_BR + 
     &                    sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BRc*Coupl23_BR )*Coupl23_A*Coupl23_Ac
     &                   /(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))

!(p1-radiation)*(p2-radiation)

	treeSS23 = treeSS23 - 2*(  (-((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad))*(  sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BLc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BRc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BLc*Coupl23_BR + 
     &                    sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BRc*Coupl23_BR )*Coupl23_A*Coupl23_Ac
     &                   /(propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2)*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))    )

!(p2-radiation)*(p2-radiation)

	treeSS23 = treeSS23 - (-(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad))*( sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BLc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BRc*Coupl23_BL - 
     &                    2*m3*m4*Coupl23_BLc*Coupl23_BR + 
     &                    sqrtS**2*(-1 - muu3**2 - muu4**2 + muu5**2 + x1 + x2)*Coupl23_BRc*Coupl23_BR)*Coupl23_A*Coupl23_Ac
     &                   /(propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                    sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))           
  

!**********************************************************Dipoles( Higgs*Higgs )**************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS23 = treeSS23 - 

!p2-splitting

!        treeSS23 = treeSS23 -  

!      endif
  
!**********************************************************************************************************************************

         enddo  
      enddo

!****************************************************Vector-exchange(gluon-exchange = 0 for elw. final states)!!!******************

      do i=sumstartV ,sumendV

        Coupl23_A = VSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
        Coupl23_BL = BL(i,f4type,f3type,f4gen,f3gen)
        Coupl23_BR = BR(i,f4type,f3type,f4gen,f3gen)
        propa = DenSmV23(i)     

         do  j=sumstartV ,sumendV 

          Coupl23_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)) * onshell_V(j) 
          Coupl23_BLc = dconjg(BL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(BR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmV23(j))

!- from Polarisation-Vectors 
       
!(p1-radiation)*(p1-radiation)
 
        treeSS23 = treeSS23 - (((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad))*( Coupl23_BL*((8*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + 
     &                    sqrtS**2*(-(muu3**4*sqrtS**2) - m4**2 - muu4**4*sqrtS**2 + muu4**2*m5**2 - 8*p1k2*x1 - 
     &                    3*m4**2*x1 + muu3**2*(8*p1k2 + sqrtS**2*(-1 + 4*muu2**2 + 2*muu4**2 + muu5**2 + x1 - 3*x2)) + 
     &                    4*m2**2*(1 + muu4**2 - muu5**2 - x1 - x2) + m4**2*x2 + 4*sqrtS**2*x1*x2))*
     &                    Coupl23_BLc + 2*m3*m4*(-4*p1k1 - 4*p1k2 + sqrtS**2*(1 - 4*muu2**2 - muu5**2 + x1 + x2))*
     &                    Coupl23_BRc) + Coupl23_BR*
     &                    (-2*m3*m4*(4*p1k1 + 4*p1k2 + sqrtS**2*(-1 + 4*muu2**2 + muu5**2 - x1 - x2))*Coupl23_BLc + 
     &                    (8*p1k1*(2*p1k2 + sqrtS**2*(muu4**2 - x2)) + sqrtS**2*
     &                    (-(muu3**4*sqrtS**2) - m4**2 - muu4**4*sqrtS**2 + muu4**2*m5**2 - 8*p1k2*x1 - 3*m4**2*x1 + 
     &                    muu3**2*(8*p1k2 + sqrtS**2*(-1 + 4*muu2**2 + 2*muu4**2 + muu5**2 + x1 - 3*x2)) + 
     &                    4*m2**2*(1 + muu4**2 - muu5**2 - x1 - x2) + m4**2*x2 + 4*sqrtS**2*x1*x2))*
     &                    Coupl23_BRc) )*Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))    

!(p1-radiation)*(p2-radiation)
 
        treeSS23 = treeSS23 - 2*( (-((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad))*( Coupl23_BL*(2*m3*m4*(1 - 2*muu1**2 - 2*muu2**2 + muu5**2)*sqrtS**2*Coupl23_BRc + 
     &                    Coupl23_BLc*(16*p1k1*p1k2 - 2*sqrtS**4 + 2*muu1**2*sqrtS**4 + 2*muu2**2*sqrtS**4 - muu3**2*sqrtS**4 + 
     &                    2*muu1**2*muu3**2*sqrtS**4 + 2*muu2**2*muu3**2*sqrtS**4 + muu3**4*sqrtS**4 - muu4**2*sqrtS**4 + 2*muu1**2*muu4**2*sqrtS**4 + 
     &                    2*muu2**2*muu4**2*sqrtS**4 - 2*muu3**2*muu4**2*sqrtS**4 + muu4**4*sqrtS**4 + 2*muu5**2*sqrtS**4 - 2*muu1**2*muu5**2*sqrtS**4 - 
     &                    2*muu2**2*muu5**2*sqrtS**4 - muu3**2*muu5**2*sqrtS**4 - muu4**2*muu5**2*sqrtS**4 - 4*p1k2*sqrtS**2*x1 + 2*sqrtS**4*x1 - 
     &                    2*muu1**2*sqrtS**4*x1 - 2*muu2**2*sqrtS**4*x1 - muu3**2*sqrtS**4*x1 + muu4**2*sqrtS**4*x1 - 4*p1k1*sqrtS**2*x2 + 
     &                    2*sqrtS**4*x2 - 2*muu1**2*sqrtS**4*x2 - 2*muu2**2*sqrtS**4*x2 + muu3**2*sqrtS**4*x2 - muu4**2*sqrtS**4*x2 
     &                    )) - Coupl23_BR*(2*m3*m4*(-1 + 2*muu1**2 + 2*muu2**2 - muu5**2)*sqrtS**2*
     &                    Coupl23_BLc + Coupl23_BRc*
     &                    (-16*p1k1*p1k2 + 2*sqrtS**4 - 2*muu1**2*sqrtS**4 - 2*muu2**2*sqrtS**4 + muu3**2*sqrtS**4 - 2*muu1**2*muu3**2*sqrtS**4 - 
     &                    2*muu2**2*muu3**2*sqrtS**4 - muu3**4*sqrtS**4 + muu4**2*sqrtS**4 - 2*muu1**2*muu4**2*sqrtS**4 - 2*muu2**2*muu4**2*sqrtS**4 + 
     &                    2*muu3**2*muu4**2*sqrtS**4 - muu4**4*sqrtS**4 - 2*muu5**2*sqrtS**4 + 2*muu1**2*muu5**2*sqrtS**4 + 2*muu2**2*muu5**2*sqrtS**4 + 
     &                    muu3**2*muu5**2*sqrtS**4 + muu4**2*muu5**2*sqrtS**4 + 4*p1k2*sqrtS**2*x1 - 2*sqrtS**4*x1 + 2*muu1**2*sqrtS**4*x1 + 
     &                    2*muu2**2*sqrtS**4*x1 + muu3**2*sqrtS**4*x1 - muu4**2*sqrtS**4*x1 + 4*p1k1*sqrtS**2*x2 - 2*sqrtS**4*x2 + 
     &                    2*muu1**2*sqrtS**4*x2 + 2*muu2**2*sqrtS**4*x2 - muu3**2*sqrtS**4*x2 + muu4**2*sqrtS**4*x2 
     &                    )) )*Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))    )

!(p2-radiation)*(p2-radiation)
 
        treeSS23 = treeSS23 - (-(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad))*(  Coupl23_BL*((-(muu3**4*sqrtS**4) + 8*p1k1*(2*p1k2 - m4**2) + 
     &                    (4*muu1**2 - muu4**2)*sqrtS**4*(1 + muu4**2 - muu5**2 - x1 - x2) + 
     &                    muu3**2*(-8*p1k2*sqrtS**2 + sqrtS**4*(-1 + 4*muu1**2 + 2*muu4**2 + muu5**2 + x1 + x2)))*Coupl23_BLc
     &                  - 2*m3*m4*(-4*p1k1 - 4*p1k2 + sqrtS**2*(-1 + 4*muu1**2 + muu5**2 + x1 + x2))*Coupl23_BRc) + 
     &                    Coupl23_BR*(2*m3*m4*(4*p1k1 + 4*p1k2 - sqrtS**2*(-1 + 4*muu1**2 + muu5**2 + x1 + x2))*
     &                    Coupl23_BLc + (-(muu3**4*sqrtS**4) + 8*p1k1*(2*p1k2 - m4**2) + 
     &                    (4*muu1**2 - muu4**2)*sqrtS**4*(1 + muu4**2 - muu5**2 - x1 - x2) + 
     &                    muu3**2*(-8*p1k2*sqrtS**2 + sqrtS**4*(-1 + 4*muu1**2 + 2*muu4**2 + muu5**2 + x1 + x2)))*Coupl23_BRc
     &                    ) )*Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                    sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))   

!**********************************************************Dipoles( Vector*Vector )************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS23 = treeSS23 - 

!p2-splitting

!        treeSS23 = treeSS23 -  

!      endif
  
!**********************************************************************************************************************************

        enddo   
  
!********************************************************Vector-Higgs *************************************************************

        do j=sumstartH ,sumendH

          Coupl23_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl23_BLc = dconjg(DL(j,f4type,f3type,f4gen,f3gen))
          Coupl23_BRc = dconjg(DR(j,f4type,f3type,f4gen,f3gen))
          propac = dconjg(DenSmH23(j))

!- from Vector-propagator * - from Polarisation-Vectors

!(p1(V)-radiation)*(p1(H)-radiation)
 
        treeSS23 = treeSS23 + 2*( (((4*p1k1 + 4*p1k2 + (-2 + 2*muu1**2 + 2*muu2**2 + muu5**2)*sqrtS**2)*V1radc*
     &                    V1rad))*(  Coupl23_BR*(m4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*
     &                    Coupl23_BLc + m3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*
     &                    Coupl23_BRc) + Coupl23_BL*
     &                    (m3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*Coupl23_BLc + 
     &                    m4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*Coupl23_BRc) )
     &                   *Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 ))      )

!(p1(V)-radiation)*(p2(H)-radiation)
 
        treeSS23 = treeSS23 + 2*( (-((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V2radc*
     &                    V1rad))*(  Coupl23_BR*(m4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*
     &                    Coupl23_BLc + m3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*
     &                    Coupl23_BRc) + Coupl23_BL*
     &                    (m3*(4*p1k2 + sqrtS**2*(-1 - muu3**2 + muu4**2 + muu5**2 + x1 - x2))*Coupl23_BLc + 
     &                    m4*(-4*p1k1 - sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 - x1 + x2))*Coupl23_BRc) )
     &                   *Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))     )

!(p2(V)-radiation)*(p1(H)-radiation)
 
        treeSS23 = treeSS23 + 2*( (-((sqrtS**2*(2*muu1**2 + 2*muu2**2 - muu5**2 - x1 - x2)*V1radc*
     &                    V2rad))*( Coupl23_BR*(m4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                    Coupl23_BLc + m3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*
     &                    Coupl23_BRc) + Coupl23_BL*
     &                    (m3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*Coupl23_BLc + 
     &                    m4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*Coupl23_BRc) )
     &                   *Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 - sqrtS**2 - m1**2 + m2**2 + m5**2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))      )

!(p2(V)-radiation)*(p2(H)-radiation)
 
        treeSS23 = treeSS23 + 2*( (-(((4*p1k1 + 4*p1k2 - sqrtS**2*(-2 + 2*muu1**2 + 2*muu2**2 + muu5**2 + 2*x1 + 2*x2))*V2radc*
     &                    V2rad))*( Coupl23_BR*(m4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*
     &                    Coupl23_BLc + m3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*
     &                    Coupl23_BRc) + Coupl23_BL*
     &                    (m3*(4*p1k2 + sqrtS**2*(1 + muu3**2 - muu4**2 - muu5**2 - x1 - x2))*Coupl23_BLc + 
     &                    m4*(-4*p1k1 + sqrtS**2*(-1 + muu3**2 - muu4**2 + muu5**2 + x1 + x2))*Coupl23_BRc) )
     &                   *Coupl23_A*Coupl23_Ac/
     &                    (propac*propa*(2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - 
     &                    sqrtS**2*x2 )*
     &                    (2*p1k1 + 2*p1k2 + sqrtS**2 - m1**2 + m2**2 - m5**2 - sqrtS**2*x1 - sqrtS**2*x2 )))      )

!**********************************************************Dipoles( Vector*Higgs )*************************************************

!      if(m5.lt.zerocheck) then    !only if radiated particle has zero-mass-dipols become necessary

!p1-splitting

!        treeSS23 = treeSS23 - 

!p2-splitting

!        treeSS23 = treeSS23 -  

!      endif
  
!**********************************************************************************************************************************
     
         enddo
       enddo

!**********************************************************************************************************************************

        QQ2xx_M2Tree23FFbar = QQ2xx_M2Tree23FFbar + NC*colorfactor23*dreal(treeSS23)

!        print *, "Tree_Mss=", NC*colorfactor23*treeSS23
!        print *, "Tree_M2Tree=", QQ2xx_M2Tree23FFbar

!==================================================================================================================================

       end
