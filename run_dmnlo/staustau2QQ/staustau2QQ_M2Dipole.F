************************************************
	double precision function staustau2QQ_M2Dipole()

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer i,j,k,l
	double precision vel,roh,roh2,roh12,roh22,VQsymm,VQasymm1,VQasymm2,GammaQ1,GammaQ2     !,Q12,Q22
	double precision mj_,sjk
	complex*16 staustau2QQ_M2Tree
	double complex cspen


        !mit 1/epsilon<->IRdiv * endliche Terme aus 1/epsilon*(4*pi)^epsilon*Gamma(1+epsilon) (für Vorfaktor siehe Loopfunctions-Datei Zeile 1)


        if((sqrtS*muu3.lt.1.d-1).and.(sqrtS*muu4.lt.1.d-1)) then

	sjk= sqrtS**2                                                                                             !sjk=2*k1*k2

          if(EpsPole.eq.-2d0) then

*singular part

          VQsymm = 1d0

	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*2d0*VQsymm
      
          elseif(EpsPole.eq.-1d0) then

*singular part

          VQsymm = dlog(muSc**2/sjk)                                                                          !hier:Vergl. Catani-Seymour (kurz) Glg. 6.20, inkl. Entwicklung des Vorfaktors 1/Gamma bis epsilon^2

* non-singular part

          VQasymm1 = 0d0                                                                                       !siehe über Glg. 6.26
 
          VQasymm2 = 0d0                                                                                       !siehe über Glg. 6.26

          GammaQ1 = 3d0/2d0*CF                                                                                 !hier:Glg. 6.28 
 
          GammaQ2 = 3d0/2d0*CF    


	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 + 1d0/CF*GammaQ1 + 1d0/CF*GammaQ2 )  !hier:3d0*dlog(muSc**2/sjk)-Term fehlt in hep/ph 9605323 Glg.7.26->nein tut er nicht; aus Vorfaktor-Entwicklung in epsilon

          elseif(EpsPole.eq.0d0) then

*singular part

          VQsymm = IRdiv**2 + IRdiv*dlog(muSc**2/sjk) + 0.5d0*dlog(muSc**2/sjk)**2 - pi**2/6d0                     !hier:Vergl. Catani-Seymour (kurz) Glg. 6.20, inkl. Entwicklung des Vorfaktors 1/Gamma bis epsilon^2

* non-singular part

          VQasymm1 = 0d0                                                                                            !siehe über Glg. 6.26
 
          VQasymm2 = 0d0                                                                                            !siehe über Glg. 6.26

          GammaQ1 = 3d0/2d0*CF*(IRdiv)-CF/2d0                                                                       !hier:Glg. 6.28 + Phys. Lett. B 500 (2001) 149-160  Glg. 17 + 20
 
          GammaQ2 = 3d0/2d0*CF*(IRdiv)-CF/2d0                                                                       !hier:Glg. 6.28 + Phys. Lett. B 500 (2001) 149-160  Glg. 17 + 20

  	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 - pi**2 + 1d0/CF*GammaQ1 + 
     &                         1d0/CF*GammaQ2 + 3d0*dlog(muSc**2/sjk) + 10d0 )                                      !hier:3d0*dlog(muSc**2/sjk)-Term fehlt in hep/ph 9605323 Glg.7.26->nein tut er nicht; aus Vorfaktor-Entwicklung in epsilon

          else

          staustau2QQ_M2Dipole= 0d0

!          print*,'EpsPole-choice not covered by Dipoles'

          endif


        elseif((sqrtS*muu3.lt.1.d-1).and.(sqrtS*muu4.ge.1.d-1)) then

        mj_=muu4*sqrtS
        sjk=sqrtS**2-mj_**2

          if(EpsPole.eq.-2d0) then

* singular part

	  VQsymm = 1d0/2d0

	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*2d0*VQsymm                                   !Vergl. Catani-Seymour (kurz) Glg.6.16

          elseif(EpsPole.eq.-1d0) then        

* singular part

	  VQsymm = dlog(muSc**2/sjk)/2d0 + 1d0/2d0*dlog(mj_**2/sjk)

* non-singular part

!j massiv

	  VQasymm1 = 0d0             !Catani-Seymour Glg. 6.22

!k massiv

          VQasymm2 = 0d0             !Catani-Seymour Glg. 6.23


          GammaQ1 = CF
          GammaQ2 = 3d0/2d0*CF

                                                                                                                                 
	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 + 1d0/CF*GammaQ1 + 1d0/CF*GammaQ2 )         

          elseif(EpsPole.eq.0d0) then   

* singular part

	  VQsymm = (IRdiv**2 + IRdiv*dlog(muSc**2/sjk) + dlog(muSc**2/sjk)**2/2d0 - pi**2/6d0)/2d0   !aber woher kommt der  dlog(sjk)-Term->aus Entwicklung des Vorfaktors nach epsilon
     &             + (IRdiv + dlog(muSc**2/sjk))/2d0*dlog(mj_**2/sjk) - dlog(mj_**2/sjk)**2/4d0 - pi**2/12d0 
     &             - dlog(mj_**2/sjk)*dlog(sjk/sqrtS**2)/2d0 - dlog(mj_**2/sqrtS**2)*dlog(sjk/sqrtS**2)/2d0         !Vergl. Catani-Seymour (kurz) Glg. 6.20 

* non-singular part

!j massiv

	  VQasymm1 = 3d0/2d0*dlog(sjk/sqrtS**2) + pi**2/6d0 - dreal(cspen(dcmplx(sjk/sqrtS**2))) - 2*dlog(sjk/sqrtS**2) - mj_**2/sjk*dlog(mj_**2/sqrtS**2)          !Catani-Seymour Glg. 6.22

!k massiv

          VQasymm2 = 3d0/2d0*(dlog(sjk/sqrtS**2) - 2*dlog((sqrtS - mj_)/sqrtS) - 2*mj_/(sqrtS + mj_) ) + pi**2/6d0 - dreal(cspen(dcmplx(sjk/sqrtS**2)))             !Catani-Seymour Glg. 6.23


          GammaQ1 = CF*(IRdiv - 1d0/2d0*dlog(muSc**2) + 1d0/2d0*dlog(mj_**2) - 2d0)
          GammaQ2 = 3d0/2d0*CF*(IRdiv)-CF/2d0                                                                       !hier:Glg. 6.28 + Phys. Lett. B 500 (2001) 149-160  Glg. 17 + 20

                                                                                                                                 
	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 - pi**2 + 1d0/CF*GammaQ1 + 
     &                         1d0/CF*GammaQ2 + 3d0*dlog(muSc**2/sjk) + 10d0 )                                      !Vergl. Catani-Seymour (kurz) Glg.6.16

          else

          staustau2QQ_M2Dipole= 0d0

!          print*,'EpsPole-choice not covered by Dipoles'

          endif

!        print*,'B'

        elseif((sqrtS*muu3.ge.1.d-1).and.(sqrtS*muu4.lt.1.d-1)) then

          mj_=muu3*sqrtS
          sjk=sqrtS**2-mj_**2

          if(EpsPole.eq.-2d0) then

* singular part

	  VQsymm = 1d0/2d0

	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*2d0*VQsymm                                   !Vergl. Catani-Seymour (kurz) Glg.6.16

          elseif(EpsPole.eq.-1d0) then        

* singular part

	  VQsymm = dlog(muSc**2/sjk)/2d0 + 1d0/2d0*dlog(mj_**2/sjk)

* non-singular part

!j massiv

	  VQasymm1 = 0d0         !Catani-Seymour Glg. 6.22

!k massiv

          VQasymm2 = 0d0             !Catani-Seymour Glg. 6.23


          GammaQ1 = CF
          GammaQ2 = 3d0/2d0*CF

                                                                                                                                 
	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 + 1d0/CF*GammaQ1 + 1d0/CF*GammaQ2 )         

          elseif(EpsPole.eq.0d0) then   

* singular part

	  VQsymm = (IRdiv**2 + IRdiv*dlog(muSc**2/sjk) + dlog(muSc**2/sjk)**2/2d0 - pi**2/6d0)/2d0   !aber woher kommt der  dlog(sjk)-Term->aus Entwicklung des Vorfaktors nach epsilon
     &             + (IRdiv + dlog(muSc**2/sjk))/2d0*dlog(mj_**2/sjk) - dlog(mj_**2/sjk)**2/4d0 - pi**2/12d0 
     &             - dlog(mj_**2/sjk)*dlog(sjk/sqrtS**2)/2d0 - dlog(mj_**2/sqrtS**2)*dlog(sjk/sqrtS**2)/2d0         !Vergl. Catani-Seymour (kurz) Glg. 6.20 

* non-singular part

!j massiv

	  VQasymm1 = 3d0/2d0*dlog(sjk/sqrtS**2) + pi**2/6d0 - dreal(cspen(dcmplx(sjk/sqrtS**2))) - 2*dlog(sjk/sqrtS**2) - mj_**2/sjk*dlog(mj_**2/sqrtS**2)          !Catani-Seymour Glg. 6.22

!k massiv

          VQasymm2 = 3d0/2d0*(dlog(sjk/sqrtS**2) - 2*dlog((sqrtS - mj_)/sqrtS) - 2*mj_/(sqrtS + mj_) ) + pi**2/6d0 - dreal(cspen(dcmplx(sjk/sqrtS**2)))             !Catani-Seymour Glg. 6.23


          GammaQ1 = CF*(IRdiv - 1d0/2d0*dlog(muSc**2) + 1d0/2d0*dlog(mj_**2) - 2d0)
          GammaQ2 = 3d0/2d0*CF*(IRdiv)-CF/2d0                                                                       !hier:Glg. 6.28 + Phys. Lett. B 500 (2001) 149-160  Glg. 17 + 20

                                                                                                                                 
	  staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 - pi**2 + 1d0/CF*GammaQ1 + 
     &                         1d0/CF*GammaQ2 + 3d0*dlog(muSc**2/sjk) + 10d0 )                                      !Vergl. Catani-Seymour (kurz) Glg.6.16

          else

          staustau2QQ_M2Dipole= 0d0

!          print*,'EpsPole-choice not covered by Dipoles'

          endif

!        print*,'C'

        elseif((sqrtS*muu3.ge.1.d-1).and.(sqrtS*muu4.ge.1.d-1)) then

          if(dabs((muu3*sqrtS-muu4*sqrtS)/muu4*sqrtS).lt.1d-6) then  !special case build in, because the more general case may become numerically unstable in the case of equal and light final-state-particles (e.g. bB)

            vel = dsqrt(1d0-4*muu3**2)
	    sjk= sqrtS**2*(1d0-2d0*muu3**2)

            if(EpsPole.eq.-2d0) then

 	    staustau2QQ_M2Dipole= 0d0

            elseif(EpsPole.eq.-1d0) then

* singular part
	    VQsymm = (1d0+vel**2)/(2d0*vel)*dlog((1d0-vel)/(1d0+vel))

	    GammaQ1 = CF

	    staustau2QQ_M2Dipole=staustau2QQ_M2Tree()*2d0*CF*GS**2/(8d0*pi**2)*(VQsymm + 1d0/CF*GammaQ1)
            
            

            elseif(EpsPole.eq.0d0) then

* singular part

	    VQsymm = (1d0+vel**2)/(2d0*vel)*(dlog((1d0-vel)/(1d0+vel))*((IRdiv + dlog(muSc**2)) - dlog(sjk)) - 
     &                         1d0/2d0*dlog((1d0-vel)/(1d0+vel))**2 - pi**2/6d0 + dlog((1d0-vel)/(1d0+vel))*dlog(2d0/(1d0+vel**2)) )

* non-singular part

	    VQasymm = - 3d0/2d0*dlog(2d0/(1d0+vel**2)) + (1d0+vel**2)/(2d0*vel)*(2d0*dlog((1d0-vel)/(1d0+vel))*dlog(2d0*(1d0+vel**2)/(1d0+vel)**2) +
     &            2d0*dreal(cspen(dcmplx((1d0-vel)**2/(1d0+vel)**2))) - 2d0*dreal(cspen(dcmplx(2d0*vel/(1d0+vel)))) -pi**2/6d0 ) + 
     &            dlog(1d0 - 1d0/2d0*dsqrt(1d0 - vel**2)) - 2d0*dlog(1d0 - dsqrt(1d0 - vel**2)) - 
     &            (1d0-vel**2)/(1d0+vel**2)*dlog(dsqrt(1d0 - vel**2)/(2d0 - dsqrt(1d0 - vel**2))) - dsqrt(1d0 - vel**2)/(2d0 - dsqrt(1d0 - vel**2)) + 
     &            2d0*(1d0 - vel**2 - dsqrt(1d0 - vel**2))/(1d0 + vel**2) + pi**2/2d0

	    GammaQ1 = CF*((IRdiv - 1d0/2d0*dlog(muSc**2)) + 1d0/2d0*dlog(sqrtS**2*muu3**2) - 2d0)

	    staustau2QQ_M2Dipole=staustau2QQ_M2Tree()*2d0*CF*GS**2/(8d0*pi**2)*(VQsymm + VQasymm - pi**2/3d0 + 1d0/CF*GammaQ1 + 3d0/2d0*dlog(muSc**2/sjk) + 5d0 - pi**2/6d0)

            else

            staustau2QQ_M2Dipole= 0d0

!            print*,'EpsPole-choice not covered by Dipoles'

            endif

          else                             

   	    sjk= sqrtS**2*(1d0-muu3**2-muu4**2)                                                                       !sjk=2*k_j*k_k
	    vel = dsqrt(1d0 - 2d0*muu3**2 - 2d0*muu4**2 + (muu3**2-muu4**2)**2)/(1d0 - muu3**2 - muu4**2)             !entspricht dem vijktilde in Catani-Seymour (kurz) Glg.5.8
            roh= dsqrt((1d0-vel)/(1d0+vel))                                                                           !Catani-Seymour (kurz) Glg.5.30
            roh2= (1d0-vel)/(1d0+vel)
            roh12= (1d0-vel+2d0*muu3**2/(1d0 - muu3**2 - muu4**2))/(1d0+vel+2d0*muu3**2/(1d0 - muu3**2 - muu4**2))    !  "              
            roh22= (1d0-vel+2d0*muu4**2/(1d0 - muu3**2 - muu4**2))/(1d0+vel+2d0*muu4**2/(1d0 - muu3**2 - muu4**2))    !  "    

            if(EpsPole.eq.-2d0) then

 	    staustau2QQ_M2Dipole= 0d0

            elseif(EpsPole.eq.-1d0) then

* singular part

	    VQsymm = 1d0/vel*dlog(roh)                                                                                !Vergl. Catani-Seymour (kurz) Glg. 6.20 

* non-singular part

	    VQasymm1 = 0d0
            VQasymm2= 0d0

	    GammaQ1 = CF                                                                                              !Catani-Seymour (kurz) Glg. 6.29.: dlog(muu3) ergibt Pol-Beitrag zu 1/epsilon für muu3->0
            GammaQ2 = CF                                                                                              !Catani-Seymour (kurz) Glg. 6.29.: dlog(muu4) ergibt Pol-Beitrag zu 1/epsilon für muu4->0
                                                                                                                                 
	    staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 + 1d0/CF*GammaQ1 +  1d0/CF*GammaQ2  )                                      !Vergl. Catani-Seymour (kurz) Glg.6.16

            elseif(EpsPole.eq.0d0) then

* singular part

	    VQsymm = 1d0/vel*(dlog(roh)*(IRdiv + dlog(muSc**2) - dlog(sjk))                                           !aber woher kommt der  dlog(sjk)-Term->aus Entwicklung des Vorfaktors nach epsilon
     &               - 1d0/4d0*dlog(roh12)**2 - 1d0/4d0*dlog(roh22)**2 - pi**2/6d0 + dlog(roh)*dlog(sqrtS**2/sjk) )   !Vergl. Catani-Seymour (kurz) Glg. 6.20 

* non-singular part

	    VQasymm1 = - 3d0/2d0*dlog(sqrtS**2/sjk) + 1d0/vel*( dlog(roh2)*dlog(1d0+roh2) +                           !Vergl. Catani-Seymour (kurz) Glg. 6.21.
     &              2d0*dreal(cspen(dcmplx(roh2))) - dreal(cspen(dcmplx(1d0-roh12)))-dreal(cspen(dcmplx(1d0-roh22))) - pi**2/6d0   )  +                  
     &              dlog(1d0-muu4) - 2d0*dlog((1d0-muu4)**2-muu3**2) - 2d0*sqrtS**2*muu3**2/sjk*dlog(muu3/(1d0-muu4))
     &              - muu4/(1d0-muu4) + 2d0*sqrtS**2*muu4*(2d0*muu4-1d0)/sjk + pi**2/2d0

            VQasymm2= - 3d0/2d0*dlog(sqrtS**2/sjk) + 1d0/vel*( dlog(roh2)*dlog(1d0+roh2) +                            !Vergl. Catani-Seymour (kurz) Glg. 6.21.
     &              2d0*dreal(cspen(dcmplx(roh2))) - dreal(cspen(dcmplx(1d0-roh12)))-dreal(cspen(dcmplx(1d0-roh22))) - pi**2/6d0   )  +                  
     &              dlog(1d0-muu3) - 2d0*dlog((1d0-muu3)**2-muu4**2) - 2d0*sqrtS**2*muu4**2/sjk*dlog(muu4/(1d0-muu3))
     &              - muu3/(1d0-muu3) + 2d0*sqrtS**2*muu3*(2d0*muu3-1d0)/sjk + pi**2/2d0

	    GammaQ1 = CF*(IRdiv - 1d0/2d0*dlog(muSc**2) + 1d0/2d0*dlog(sqrtS**2*muu3**2) - 2d0)                       !Catani-Seymour (kurz) Glg. 6.29.: dlog(muu3) ergibt Pol-Beitrag zu 1/epsilon für muu3->0
            GammaQ2 = CF*(IRdiv - 1d0/2d0*dlog(muSc**2) + 1d0/2d0*dlog(sqrtS**2*muu4**2) - 2d0)                       !Catani-Seymour (kurz) Glg. 6.29.: dlog(muu4) ergibt Pol-Beitrag zu 1/epsilon für muu4->0
                                                                                                                                 
	    staustau2QQ_M2Dipole= staustau2QQ_M2Tree()*CF*GS**2/(8d0*pi**2)*(2d0*VQsymm + VQasymm1 + VQasymm2 - pi**2 + 1d0/CF*GammaQ1 + 
     &                           1d0/CF*GammaQ2 + 3d0*dlog(muSc**2/sjk) + 10d0 )                                      !Vergl. Catani-Seymour (kurz) Glg.6.16

            else

            staustau2QQ_M2Dipole= 0d0



            endif

          endif


                                                                                                                                  
        else

        staustau2QQ_M2Dipole= 0d0

        print*,'case not covered by integrated Dipoles'

        endif

                                                                                                                                    

	end
