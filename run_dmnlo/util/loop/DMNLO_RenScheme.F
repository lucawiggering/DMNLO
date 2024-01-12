        integer function check_RenScheme(renscheme, renscale, slhafilename)
        implicit none
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_Model.h"
        integer MOswitch,renscheme
        character slhafilename*200
        double precision renscale

        QScale = renscale
        Qscalealphas = Qscale

        MOswitch = 0
        call DMNLO_ModelPara(MOswitch,slhafilename)
        call DMNLO_ModelIni(slhafilename)

        choosesol = 0 ! the check is only contained for this option so far 
        call Init_RenScheme(renscheme)

        check_RenScheme = renfail

      end function 
      
      subroutine Init_RenScheme(renscheme)
      implicit none

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

      integer i,j,k,l,m,n,type,gen, renscheme
      integer idv1,idv2,idmv1,idmv2,iv1,iv2
      integer nf
      
      double precision M2LL,M2RR,M2LR

      double precision eb,et
      double precision dtht(4,3),dthtm1m2(4,3),dM2LR,dmfdr,dmsf12,dmsf22
      double precision topLL, topLR, topRR, botLL, botLR, botRR
      double precision MQ2v1,MD2v1,MU2v1,MQ2v2,MD2v2,MU2v2
      double precision MSf2temp(2,4,3)
      
      double precision u11,u12,u21,u22,ut11,ut12,ut21,ut22,ub11,ub12,ub21,ub22
      double precision mstv1(2),mstv12(2),mstv2(2),mstv22(2),msb(2),msb2(2)
      double precision stv1,stv2,sttv1,sttv2,sbv1,sbv2
      double complex ustv1(2,2),ustv2(2,2),usbv1(2,2),usbv2(2,2)

      
* functions      
      double precision ASDR,MBDR,MTDR,Epsij,MB2OS
      
* subroutines
      external DMNLO_CalcPif, DMNLO_CalcdPif, DMNLO_CalcPiSf, DMNLO_CalcdPiSf, DMNLO_CalcdPiG

! renfail is a flag that determines whether the renormalisation scheme works or not (1 - successful, 0 - not successful)
      renfail = 1

! ===== scale dependencies =====	
!     QScale should be read out of SPheno/SoftSUSY to be consistent with all our DRbar parameters


      muSc = QScale


! ===== Running alpha(s) in DR & MSSM =====

      Alfas = (4d0*pi)*ASDR(Qscale**2)
      Alfas2 = Alfas**2
      GS = sqrt(4*pi*Alfas)

! ===== Top and Bottom charges =====

      et =  2./3.
      eb = -1./3.

! ===== initialie dUSF =====

      do type=3,4
      do gen=1,3

      dUSf(1,1,type,gen) = (0d0,0d0)
      dUSf(1,2,type,gen) = (0d0,0d0)
      dUSf(2,1,type,gen) = (0d0,0d0)
      dUSf(2,2,type,gen) = (0d0,0d0)

      enddo
      enddo

! =====================================
! =========== DR Scheme Input =========
! ==                                 ==
! == mb,mt,MSf,theta,Af all in DRbar ==
! ======================== 

      if (renscheme.eq.0) then

* Quark DRbar masses
      Mf(3,3) = MTDR(MTos,Qscale)
      Mf(4,3) = MBDR(MBSM,MBSM**2,Qscale**2,5,3)

      Mf2(3,3) = Mf(3,3)**2
      Mf2(4,3) = Mf(4,3)**2
		 
      do gen=1,3 
		
* up sector DRbar         
      M2LL = MSQ(gen)**2 + Mf(3,gen)**2 + C2B*mz**2*(0.5 - et*SW2)
      M2RR = MSU(gen)**2 + Mf(3,gen)**2 + C2B*et*mz**2*SW2
      M2LR = Mf(3,gen)*(Af(3,gen) - dreal(MUE)/TB)
      call Eigen2x2(MSf(1,3,gen), MSf2(1,3,gen), USf(1,1,3,gen), M2LL,M2RR,dcmplx(M2LR), *999)


* down sector DRbar         
      M2LL = MSQ(gen)**2 + Mf(4,gen)**2 + C2B*mz**2*(-0.5 - eb*SW2)
      M2RR = MSD(gen)**2 + Mf(4,gen)**2 + C2B*eb*mz**2*SW2
      M2LR = Mf(4,gen)*(Af(4,gen) - dreal(MUE)*TB)
      call Eigen2x2(MSf(1,4,gen), MSf2(1,4,gen), USf(1,1,4,gen), M2LL,M2RR,dcmplx(M2LR), *999)
			
				 
      enddo
		
	
! =====================================
! ======= DR Scheme Counterterms ======
! ===================================== 

*************************************************
* Alpha_s renormalization constant (DRbar, MSSM)
* deltaGs is in fact dGS/GS
         nf = 6
         DeltaGs = 0d0
         
! usual version         
         DeltaGs = (Alfas/(8d0*pi))*(UVdiv*(nf - 3*CV))
                  
! debugging version - splitting into SM and SUSY contrib

! ! SM contrib - gluon + ghost + quark loop
         DeltaGs = 0d0
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(2d0*nf/3d0 - 11d0*CV/3d0))
! ! SUSY contrib - gluino loop
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(1d0*nf/3d0 + 2d0*CV/3d0))

      
*************************************************
* Gluon wave-function renormalization constant

      call DMNLO_CalcdPiG(0d0)
      dZG = - dreal(dPiGT)

! Ghost wave-function renormalization constant
      call DMNLO_CalcdPiC(0d0)
      dZC=-dreal(dPiC)
   


***************************
* Quark mass counterterms
        do type=3,4 
          do gen=1,3
                dMf(type,gen) = Mf(type,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)
          enddo    
        enddo          

*************************************************
* Quark wave-function renormalization constants
	  do type=3,4
	    do gen=1,3
	        call DMNLO_CalcPif(Mf2(type,gen),type,gen)
	        call DMNLO_CalcdPif(Mf2(type,gen),type,gen)
	        dZfL(type,gen) = -PifVL(type,gen) - Mf2(type,gen)*(dPifVL(type,gen) + dPifVR(type,gen)) - Mf(type,gen)*dPifSL(type,gen) - Mf(type,gen)*dPifSR(type,gen)
	        dZfR(type,gen) = -PifVR(type,gen) - Mf2(type,gen)*(dPifVR(type,gen) + dPifVL(type,gen)) - Mf(type,gen)*dPifSR(type,gen) - Mf(type,gen)*dPifSL(type,gen)
	    enddo
	  enddo


***************************
* Squark mass counterterms
        do type=3,4
          do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))
         
              dMSf2(1,type,gen) = GS**2*CF/(16d0*pi**2)*( (u11**2-u12**2)**2*MSf2(1,type,gen) + (u21*u11-u22*u12)**2*MSf2(2,type,gen) - 
     &                            MSf2(1,type,gen) + 8d0*(u11*u12)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv
         
              dMSf2(2,type,gen) = GS**2*CF/(16d0*pi**2)*( (u21**2-u22**2)**2*MSf2(2,type,gen) + (u21*u11-u22*u12)**2*MSf2(1,type,gen) - 
     &                            MSf2(2,type,gen) + 8d0*(u21*u22)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv
	    enddo
	  enddo

**************************************************
* Squark wave-function renormalization constants
	  do type=3,4
	    do gen=1,3
	      do k=1,2
	        do l=1,2
	          if (k.eq.l) then
	           call DMNLO_CalcdPiSf(MSf2(k,type,gen),type,gen)
	           dZSf(k, k, type,gen) = -dPiSf(k,k,type,gen)
	          else 
	           call DMNLO_CalcPiSf(MSf2(l,type,gen),type,gen)
	           dZSf(k, l, type,gen) = (2d0/(MSf2(k,type,gen) - MSf2(l,type,gen)))*PiSf(k,l,type,gen) 
	          endif
	        enddo
	      enddo
	    enddo
	  enddo
        
! **************************************************
!  Squark mixing angle & mixing matrix counterterms

        do type=3,4
          do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))
        
              dtht(type,gen) = GS**2*CF/(16d0*pi**2)/(MSf2(1,type,gen)- MSf2(2,type,gen))*
     &               ( (u11*u21-u22*u12)*((u11**2-u12**2)*MSf2(1,type,gen) + (u21**2-u22**2)*MSf2(2,type,gen)) + 
     &               4d0*(u11*u22+u12*u21)*MGl*Mf(type,gen))*UVdiv
      
              do m=1,2
                dUSf(1,1,type,gen) = dUSf(1,1,type,gen) + Epsij(1,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(1,2,type,gen) = dUSf(1,2,type,gen) + Epsij(1,m)*USf(m,2,type,gen)*dtht(type,gen)
                dUSf(2,1,type,gen) = dUSf(2,1,type,gen) + Epsij(2,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(2,2,type,gen) = dUSf(2,2,type,gen) + Epsij(2,m)*USf(m,2,type,gen)*dtht(type,gen)
              enddo
              
	    enddo
	  enddo
! ***************************
! Af counterterms

        do type=3,4
          do gen=1,3
             u11=dreal(USf(1,1,type,gen))
             u12=dreal(USf(1,2,type,gen))
             u21=dreal(USf(2,1,type,gen))
             u22=dreal(USf(2,2,type,gen))
        
             dthtm1m2(type,gen) = GS**2*CF/(16d0*pi**2)*
     &              ( (u11*u21-u22*u12)*((u11**2-u12**2)*MSf2(1,type,gen) + (u21**2-u22**2)*MSf2(2,type,gen)) + 
     &               4d0*(u11*u22+u12*u21)*MGl*Mf(type,gen))*UVdiv
        
      
             dM2LR = u11*u12*dMSf2(1,type,gen) + u21*u22*dMSf2(2,type,gen)
     &               + (u21*u12 + u11*u22)*dthtm1m2(type,gen)
             
             dAf(type,gen) = 1d0/Mf(type,gen)*(dM2LR - dMf(type,gen)/Mf(type,gen)*(u11*u12*MSf2(1,type,gen) + u21*u22*MSf2(2,type,gen)))
        
          enddo
        enddo
      
! =====================================================
! ======= Mixed DR-OS Scheme No.1 - Input Param. ======
! =======                                        ======
! ============ mb,Ab,At defined in DRbar ==============
! =========  ========
! ========= theta, mst2 dependant parameters ==========
! ===================================================== 

      else if (renscheme.eq.1) then
      


***************************
* Bottom quark DRbar mass
            Mf(4,3) = MBDR(MBSM,MBSM**2,Qscale**2,5,3)
            Mf2(4,3) = Mf(4,3)**2
           
***************************
* Top quark pole mass
            Mf(3,3) = Mfos(3,3)
            Mf2(3,3) = Mf(3,3)**2

***************************************
* Stop & Sbottom sector from msb1,msb2

* solution 1 
           MQ2v1 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 + C2B*MZ**2 + 2d0*C2B*eb*MZ**2*SW2 - 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           MD2v1 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 - 2d0*C2B*eb*MZ**2*SW2 + 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           
           MU2v1 = (-2d0*MQ2v1*MSf(1,3,3)**2 + 2*MSf(1,3,3)**4 - 2*Af(3,3)**2*Mf(3,3)**2 + 2*MQ2v1*Mf(3,3)**2 - 
     &         4*MSf(1,3,3)**2*Mf(3,3)**2 + 2*Mf(3,3)**4 - 
     &         c2b*MSf(1,3,3)**2*mz**2 + c2b*Mf(3,3)**2*mz**2 + 2*c2b*et*MQ2v1*mz**2*SW2 + c2b**2*et*mz**4*SW2 - 
     &         2*c2b**2*et**2*mz**4*SW2**2 - (2*Mf(3,3)**2*MUE**2)/tb**2+ (4*Af(3,3)*Mf(3,3)**2*MUE)/tb)/
     &       (-2*MQ2v1 + 2*MSf(1,3,3)**2 - 2*Mf(3,3)**2 - c2b*mz**2 + 2*c2b*et*mz**2*SW2)
           
           topLL = MQ2v1 + Mf(3,3)**2 + C2B*MZ**2*(0.5d0 - et*SW2)
           topRR = MU2v1 + Mf(3,3)**2 + C2B*et*MZ**2*SW2
           topLR = Mf(3,3)*(Af(3,3) - MUE/TB)
                      
           call Eigen2x2(mstv1(1), mstv12(1), ustv1(1,1),
     &             topLL,topRR,dcmplx(topLR), *999)

           botLL = MQ2v1 + Mf(4,3)**2 + C2B*MZ**2*(-0.5d0 - eb*SW2)
           botRR = MD2v1 + Mf(4,3)**2 + C2B*eb*MZ**2*SW2
           botLR = Mf(4,3)*(Af(4,3) - MUE*TB)

           call Eigen2x2(msb(1), msb2(1), usbv1(1,1),
     &             botLL,botRR,dcmplx(botLR), *999)
            	

! solution 2 
           MQ2v2 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 + C2B*MZ**2 + 2d0*C2B*eb*MZ**2*SW2 + 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           MD2v2 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 - 2d0*C2B*eb*MZ**2*SW2 - 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           
           MU2v2 = (-2d0*MQ2v2*MSf(1,3,3)**2 + 2*MSf(1,3,3)**4 - 2*Af(3,3)**2*Mf(3,3)**2 + 2*MQ2v2*Mf(3,3)**2 - 
     &         4*MSf(1,3,3)**2*Mf(3,3)**2 + 2*Mf(3,3)**4 - 
     &         c2b*MSf(1,3,3)**2*mz**2 + c2b*Mf(3,3)**2*mz**2 + 2*c2b*et*MQ2v2*mz**2*SW2 + c2b**2*et*mz**4*SW2 - 
     &         2*c2b**2*et**2*mz**4*SW2**2 - (2*Mf(3,3)**2*MUE**2)/tb**2+ (4*Af(3,3)*Mf(3,3)**2*MUE)/tb)/
     &       (-2*MQ2v2 + 2*MSf(1,3,3)**2 - 2*Mf(3,3)**2 - c2b*mz**2 + 2*c2b*et*mz**2*SW2)
           
           topLL = MQ2v2 + Mf(3,3)**2 + C2B*MZ**2*(0.5d0 - et*SW2)
           topRR = MU2v2 + Mf(3,3)**2 + C2B*et*MZ**2*SW2
           topLR = Mf(3,3)*(Af(3,3) - MUE/TB)
           
           call Eigen2x2(mstv2(1), mstv22(1), ustv2(1,1),
     &             topLL,topRR,dcmplx(topLR), *999)

           botLL = MQ2v1 + Mf(4,3)**2 + C2B*MZ**2*(-0.5d0 - eb*SW2)
           botRR = MD2v1 + Mf(4,3)**2 + C2B*eb*MZ**2*SW2
           botLR = Mf(4,3)*(Af(4,3) - MUE*TB)

           call Eigen2x2(msb(1), msb2(1), usbv2(1,1),
     &             botLL,botRR,dcmplx(botLR), *999)



            if (choosesol.eq.1) then !Program chooses which solution should be taken
                           do m=1,2
                                 MSf(m,3,3) = mstv1(m)
                                 MSf2(m,3,3) = mstv12(m)
                                 do n=1,2
                                   USf(m,n,3,3) = ustv1(m,n)
                                   USf(m,n,4,3) = usbv1(m,n)
                                 enddo  
                           enddo              
            elseif (choosesol.eq.2) then !Program chooses which solution should be taken
                        do m=1,2
                              MSf(m,3,3) = mstv2(m)
                              MSf2(m,3,3) = mstv22(m)
                              do n=1,2
                                USf(m,n,3,3) = ustv2(m,n)
                                USf(m,n,4,3) = usbv2(m,n)
                              enddo  
                        enddo              
            elseif (choosesol.eq.0) then !Program chooses which solution should be taken 

* testing if any of the solutions satisfies all conditions to define counterterms in a stable manner
* 1. condition - dependant mass has to be the heavier one
* 2. condition - mixing for both stops & sbottoms 
!original version

!             idmv1 = 1
!             idmv2 = 1
!             idv1 = 1
!             idv2 = 1
!             
!             if ((mstv1(1)-MSf(1,3,3))/MSf(1,3,3).gt.1.d-5) idmv1=0
!             if ((mstv2(1)-MSf(1,3,3))/MSf(1,3,3).gt.1.d-5) idmv2=0

!new version, added by Moritz to sort out potential NANs, which can appear e.g. inMQv1,MQv2 for special parameter-combinations (->negative argument in sqrt)

            idmv1 = 0
            idmv2 = 0
            idv1 = 1
            idv2 = 1
            
            if (dabs(mstv1(1)-MSf(1,3,3))/MSf(1,3,3).le.1.d-5) idmv1=1   !MO:should sort out potential NANs
            if (dabs(mstv2(1)-MSf(1,3,3))/MSf(1,3,3).le.1.d-5) idmv2=1
            
            stv1 = dabs(dreal(ustv1(1,1)*ustv1(1,2)))
            sbv1 = dabs(dreal(usbv1(1,1)*usbv1(1,2)))
            stv2 = dabs(dreal(ustv2(1,1)*ustv2(1,2)))
            sbv2 = dabs(dreal(usbv2(1,1)*usbv2(1,2)))

            sttv1 = dabs(dreal(ustv1(1,1)*ustv1(2,2) + ustv1(1,2)*ustv1(2,1)))
            sttv2 = dabs(dreal(ustv2(1,1)*ustv2(2,2) + ustv2(1,2)*ustv2(2,1)))

            if ((stv1.le.1d-4).or.(sttv1.le.1d-4).or.(sbv1.le.1d-4)) idv1 = 0
            if ((stv2.le.1d-4).or.(sttv2.le.1d-4).or.(sbv2.le.1d-4)) idv2 = 0

            iv1 = idv1*idmv1
            iv2 = idv2*idmv2

            if ((iv1.eq.0).and.(iv2.eq.0)) then
                    print*,'*** Warning: Both solutions incompatible ***'
                    renfail = 0
            else
*           solution 1 incompatible
                  if (iv1.eq.0) then
                        do m=1,2
                              MSf(m,3,3) = mstv2(m)
                              MSf2(m,3,3) = mstv22(m)
                              do n=1,2
                                USf(m,n,3,3) = ustv2(m,n)
                                USf(m,n,4,3) = usbv2(m,n)
                              enddo  
                        enddo      
*           solution 2 incompatible
                  else if (iv2.eq.0) then
                        do m=1,2
                              MSf(m,3,3) = mstv1(m)
                              MSf2(m,3,3) = mstv12(m)
                              do n=1,2
                                USf(m,n,3,3) = ustv1(m,n)
                                USf(m,n,4,3) = usbv1(m,n)
                              enddo  
                        enddo                                    
*           both solutions compatible
                  else
                     if ( dabs(mstv1(2)-MSf(2,3,3)).lt.dabs(mstv2(2)-MSf(2,3,3)) ) then
                           do m=1,2
                                 MSf(m,3,3) = mstv1(m)
                                 MSf2(m,3,3) = mstv12(m)
                                 do n=1,2
                                   USf(m,n,3,3) = ustv1(m,n)
                                   USf(m,n,4,3) = usbv1(m,n)
                                 enddo  
                           enddo                                       
                     else
                           do m=1,2
                                 MSf(m,3,3) = mstv2(m)
                                 MSf2(m,3,3) = mstv22(m)
                                 do n=1,2
                                   USf(m,n,3,3) = ustv2(m,n)
                                   USf(m,n,4,3) = usbv2(m,n)
                                 enddo  
                           enddo                                       
                     endif      
                  endif
            endif
            
            endif !End of choosesol
            USf(2,1,3,3) = USf(2,1,3,3)*(-1)
            USf(2,2,3,3) = USf(2,2,3,3)*(-1)
            
            



! =====================================================
! ======= Mixed DR-OS Scheme No.1 - Counterterms ======
! =======                                        ======
! ============ mb,Ab,At defined in DRbar ==============
! ========= mt,mst1,msb1,msb2 defined on-shell ========
! ========= theta, mst2 dependant parameters ==========
! ===================================================== 

*************************************************
* Alpha_s renormalization constant (DRbar, MSSM)
* deltaGs is in fact dGS/GS
         nf = 6
         DeltaGs = 0d0
         
! usual version         
         DeltaGs = (Alfas/(8d0*pi))*(UVdiv*(nf - 3*CV))
                  
! debugging version - splitting into SM and SUSY contrib

! ! SM contrib - gluon + ghost + quark loop
         DeltaGs = 0d0
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(2d0*nf/3d0 - 11d0*CV/3d0))
! ! SUSY contrib - gluino loop
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(1d0*nf/3d0 + 2d0*CV/3d0))


     
*************************************************
* Gluon wave-function renormalization constant

      call DMNLO_CalcdPiG(0d0)	        
      dZG = - dreal(dPiGT)     
   
! Ghost wave-function renormalization constant
      call DMNLO_CalcdPiC(0d0)
      dZC=-dreal(dPiC)
	 
	  
***************************
* Gluino mass counterterms
      call DMNLO_CalcPigl(MGl**2)
      dMGl = MGl/(2d0)*(PiglVL + PiglVR + PiglSL/mGl + PiglSR/mGl)


*************************************************
* Gluino wave-function renormalization constants

      call DMNLO_CalcPigl(MGl**2)
      call DMNLO_CalcdPigl(MGl**2)
      dZGlL = -PiglVL - MGl**2*(dPiglVL + dPiglVR + dPiglSL/mGl + dPiglSR/mGl)
      dZGlR = -PiglVR - MGl**2*(dPiglVR + dPiglVL + dPiglSR/mGl + dPiglSL/mGl)
	        
     
**************************            
* Quark mass counterterms

* DR counterterms for all quarks but the top quark
        do type=3,4 
          do gen=1,3
             dMf(type,gen) = Mf(type,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)
          enddo    
        enddo          
         
* OS counterterm for the top quark
        call DMNLO_CalcPif(Mf2(3,3),3,3)
	  dMf(3,3) = (Mf(3,3)/2d0)*(PifVL(3,3) + PifVR(3,3)) + (PifSL(3,3) + PifSR(3,3))/2d0

***********************************************
* Quark wave-function renormalization constants
	  do type=3,4
	    do gen=1,3
	        call DMNLO_CalcPif(Mf2(type,gen),type,gen)
	        call DMNLO_CalcdPif(Mf2(type,gen),type,gen)
	        dZfL(type,gen) = -PifVL(type,gen) - Mf2(type,gen)*(dPifVL(type,gen) + dPifVR(type,gen)) - Mf(type,gen)*dPifSL(type,gen) - Mf(type,gen)*dPifSR(type,gen)
	        dZfR(type,gen) = -PifVR(type,gen) - Mf2(type,gen)*(dPifVR(type,gen) + dPifVL(type,gen)) - Mf(type,gen)*dPifSR(type,gen) - Mf(type,gen)*dPifSL(type,gen)
	    enddo
	  enddo

***************************
* Squark mass counterterms

* OS for sbottom 1,2 and stop 1
	do type=3,4
	  do gen=1,3
	    do k=1,2
	
	      call DMNLO_CalcPiSf(MSf2(k,type,gen),type,gen)
        	dMSf2(k,type,gen) = PiSf(k,k,type,gen)
	
	    enddo
	  enddo
	enddo

* stop 2 dependant counterterm - defined BELOW

**************************************************
* Squark wave-function renormalization constants

	  do type=3,4
	    do gen=1,3
	      do k=1,2
	        do l=1,2
	          if (k.eq.l) then
	           call DMNLO_CalcdPiSf(MSf2(k,type,gen),type,gen)
	           dZSf(k, k, type,gen) = -dPiSf(k,k,type,gen)
	          else 
	           call DMNLO_CalcPiSf(MSf2(l,type,gen),type,gen)
	           dZSf(k, l, type,gen) = (2d0/(MSf2(k,type,gen) - MSf2(l,type,gen)))*PiSf(k,l,type,gen) 
	          endif
	        enddo
	      enddo
	    enddo
	  enddo

***************************
* Af counterterms in DRbar
      
         do type=3,4
           do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))

              dmfdr = Mf(type,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)            

              dthtm1m2(type,gen) = GS**2*CF/(16d0*pi**2)*
     &              ( (u11*u21-u22*u12)*((u11**2-u12**2)*MSf2(1,type,gen) + (u21**2-u22**2)*MSf2(2,type,gen)) + 
     &               4d0*(u11*u22+u12*u21)*MGl*Mf(type,gen))*UVdiv

              dmsf12 = GS**2*CF/(16d0*pi**2)*( (u11**2-u12**2)**2*MSf2(1,type,gen) + (u21*u11-u22*u12)**2*MSf2(2,type,gen) - 
     &                            MSf2(1,type,gen) + 8d0*(u11*u12)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv
    
              dmsf22 = GS**2*CF/(16d0*pi**2)*( (u21**2-u22**2)**2*MSf2(2,type,gen) + (u21*u11-u22*u12)**2*MSf2(1,type,gen) - 
     &                            MSf2(2,type,gen) + 8d0*(u21*u22)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv


              dM2LR = u11*u12*dmsf12 + u21*u22*dmsf22 + (u21*u12 + u11*u22)*dthtm1m2(type,gen)

              dAf(type,gen) = 1d0/Mf(type,gen)*(dM2LR - dmfdr/Mf(type,gen)*(u11*u12*MSf2(1,type,gen) + u21*u22*MSf2(2,type,gen)))

           enddo
         enddo
      
***************************************************
* Squark mixing angle & mixing matrix counterterms

         do type=3,4
           do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))
              
              if (type.eq.3) then                              
                  dM2LR = dMf(type,gen)*(Af(type,gen)-MUE/TB) + Mf(type,gen)*dAf(type,gen)
              else 
                  dM2LR = dMf(type,gen)*(Af(type,gen)-MUE*TB) + Mf(type,gen)*dAf(type,gen)
              endif
              
              dtht(type,gen) = 1d0/((u21*u12 + u11*u22)*(MSf2(1,type,gen)-MSf2(2,type,gen)))*(dM2LR - u11*u12*(dMSf2(1,type,gen) - dMSf2(2,type,gen)))
              dthtm1m2(type,gen) = 1d0/(u21*u12 + u11*u22)*(dM2LR - u11*u12*(dMSf2(1,type,gen) - dMSf2(2,type,gen)))
              
               do m=1,2
                do n=1,2
                   dUSf(m,n,4,3) = 0d0
                enddo
              enddo


              do m=1,2
                dUSf(1,1,type,gen) = dUSf(1,1,type,gen) + Epsij(1,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(1,2,type,gen) = dUSf(1,2,type,gen) + Epsij(1,m)*USf(m,2,type,gen)*dtht(type,gen)
                dUSf(2,1,type,gen) = dUSf(2,1,type,gen) + Epsij(2,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(2,2,type,gen) = dUSf(2,2,type,gen) + Epsij(2,m)*USf(m,2,type,gen)*dtht(type,gen)
              enddo
				
           enddo
         enddo

***************************************************
* Dependant counterterms dtheta, dUSf, dmst2

         ut11=dreal(USf(1,1,3,3))
         ut12=dreal(USf(1,2,3,3))
         ut21=dreal(USf(2,1,3,3))
         ut22=dreal(USf(2,2,3,3))
         
         ub11=dreal(USf(1,1,4,3))
         ub12=dreal(USf(1,2,4,3))
         ub21=dreal(USf(2,1,4,3))
         ub22=dreal(USf(2,2,4,3))
         
            
         dMSf2(2,3,3) = (ut21*ut12 + ut11*ut22)/(ut21*ut12)*
     &   (ub11**2*dMSf2(1,4,3) + ub21**2*dMSf2(2,4,3) + 2d0*ub11*ub21*dthtm1m2(4,3) - 2d0*Mf(4,3)*dMf(4,3)
     &   -ut11**2*dMSf2(1,3,3) + 2d0*Mf(3,3)*dMf(3,3)) -
     &   2d0*ut11/ut12*(dMf(3,3)*(Af(3,3)-MUE/TB) + Mf(3,3)*dAf(3,3) - ut11*ut12*dMSf2(1,3,3))          

         dtht(3,3) = 1d0/((ut21*ut12 + ut11*ut22)*(MSf2(1,3,3)-MSf2(2,3,3)))*
     &               (dMf(3,3)*(Af(3,3)-MUE/TB) + Mf(3,3)*dAf(3,3) - ut11*ut12*(dMSf2(1,3,3) - dMSf2(2,3,3)))

        do m=1,2
         do n=1,2
            dUSf(m,n,3,3) = 0d0
         enddo
       enddo

         do m=1,2
           dUSf(1,1,3,3) = dUSf(1,1,3,3) + Epsij(1,m)*USf(m,1,3,3)*dtht(3,3)
           dUSf(1,2,3,3) = dUSf(1,2,3,3) + Epsij(1,m)*USf(m,2,3,3)*dtht(3,3)
           dUSf(2,1,3,3) = dUSf(2,1,3,3) + Epsij(2,m)*USf(m,1,3,3)*dtht(3,3)
           dUSf(2,2,3,3) = dUSf(2,2,3,3) + Epsij(2,m)*USf(m,2,3,3)*dtht(3,3)
         enddo

         
! =====================================================
! ======= Mixed DR-OS Scheme No.2 - Input Param. ======
! =======                                        ======
! ============ mt,mb,Ab,At defined in DRbar ===========
! ========= mst1,msb1,msb2 defined on-shell ===========
! ========= theta, mst2 dependant parameters ==========
! ===================================================== 

      else if (renscheme.eq.2) then
      


***************************
* Bottom quark DRbar mass
           Mf(4,3) = MBDR(MBSM,MBSM**2,Qscale**2,5,3)
           Mf2(4,3) = Mf(4,3)**2
           
***************************
* Top DRbar mass
           Mf(3,3) = MTDR(MTos,Qscale)
           Mf2(3,3) = Mf(3,3)**2

***************************************
* Stop & Sbottom sector from msb1,msb2

* solution 1 
           MQ2v1 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 + C2B*MZ**2 + 2d0*C2B*eb*MZ**2*SW2 - 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           MD2v1 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 - 2d0*C2B*eb*MZ**2*SW2 + 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           
           MU2v1 = (-2d0*MQ2v1*MSf(1,3,3)**2 + 2*MSf(1,3,3)**4 - 2*Af(3,3)**2*Mf(3,3)**2 + 2*MQ2v1*Mf(3,3)**2 - 
     &         4*MSf(1,3,3)**2*Mf(3,3)**2 + 2*Mf(3,3)**4 - 
     &         c2b*MSf(1,3,3)**2*mz**2 + c2b*Mf(3,3)**2*mz**2 + 2*c2b*et*MQ2v1*mz**2*SW2 + c2b**2*et*mz**4*SW2 - 
     &         2*c2b**2*et**2*mz**4*SW2**2 - (2*Mf(3,3)**2*MUE**2)/tb**2 + (4*Af(3,3)*Mf(3,3)**2*MUE)/tb)/
     &       (-2*MQ2v1 + 2*MSf(1,3,3)**2 - 2*Mf(3,3)**2 - c2b*mz**2 + 2*c2b*et*mz**2*SW2)
           
           topLL = MQ2v1 + Mf(3,3)**2 + C2B*MZ**2*(0.5d0 - et*SW2)
           topRR = MU2v1 + Mf(3,3)**2 + C2B*et*MZ**2*SW2
           topLR = Mf(3,3)*(Af(3,3) - MUE/TB)
                      
           call Eigen2x2(mstv1(1), mstv12(1), ustv1(1,1),
     &             topLL,topRR,dcmplx(topLR), *999)

           botLL = MQ2v1 + Mf(4,3)**2 + C2B*MZ**2*(-0.5d0 - eb*SW2)
           botRR = MD2v1 + Mf(4,3)**2 + C2B*eb*MZ**2*SW2
           botLR = Mf(4,3)*(Af(4,3) - MUE*TB)

           call Eigen2x2(msb(1), msb2(1), usbv1(1,1),
     &             botLL,botRR,dcmplx(botLR), *999)
            

* solution 2 
           MQ2v2 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 + C2B*MZ**2 + 2d0*C2B*eb*MZ**2*SW2 + 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           MD2v2 = (-2d0*Mf(4,3)**2 + MSf(1,4,3)**2 + MSf(2,4,3)**2 - 2d0*C2B*eb*MZ**2*SW2 - 
     &         dsqrt(dreal(MSf(1,4,3)**4 - 4d0*Af(4,3)**2*Mf(4,3)**2 - 2*MSf(1,4,3)**2*MSf(2,4,3)**2 + MSf(2,4,3)**4 + 
     &         8d0*Af(4,3)*Mf(4,3)**2*MUE*tb - 4d0*Mf(4,3)**2*MUE**2*tb**2)))/2d0
           
           
           MU2v2 = (-2d0*MQ2v2*MSf(1,3,3)**2 + 2*MSf(1,3,3)**4 - 2*Af(3,3)**2*Mf(3,3)**2 + 2*MQ2v2*Mf(3,3)**2 - 
     &         4*MSf(1,3,3)**2*Mf(3,3)**2 + 2*Mf(3,3)**4 - 
     &         c2b*MSf(1,3,3)**2*mz**2 + c2b*Mf(3,3)**2*mz**2 + 2*c2b*et*MQ2v2*mz**2*SW2 + c2b**2*et*mz**4*SW2 - 
     &         2*c2b**2*et**2*mz**4*SW2**2 - (2*Mf(3,3)**2*MUE**2)/tb**2 + (4*Af(3,3)*Mf(3,3)**2*MUE)/tb)/
     &       (-2*MQ2v2 + 2*MSf(1,3,3)**2 - 2*Mf(3,3)**2 - c2b*mz**2 + 2*c2b*et*mz**2*SW2)
           
           topLL = MQ2v2 + Mf(3,3)**2 + C2B*MZ**2*(0.5d0 - et*SW2)
           topRR = MU2v2 + Mf(3,3)**2 + C2B*et*MZ**2*SW2
           topLR = Mf(3,3)*(Af(3,3) - MUE/TB)
           
           call Eigen2x2(mstv2(1), mstv22(1), ustv2(1,1),
     &             topLL,topRR,dcmplx(topLR), *999)

           botLL = MQ2v1 + Mf(4,3)**2 + C2B*MZ**2*(-0.5d0 - eb*SW2)
           botRR = MD2v1 + Mf(4,3)**2 + C2B*eb*MZ**2*SW2
           botLR = Mf(4,3)*(Af(4,3) - MUE*TB)

           call Eigen2x2(msb(1), msb2(1), usbv2(1,1),
     &             botLL,botRR,dcmplx(botLR), *999)

           

* testing if any of the solutions satisfies all conditions to define counterterms in a stable manner
* 1. condition - dependant mass has to be the heavier one
* 2. condition - mixing for both stops & sbottoms 
!original version

!             idmv1 = 1
!             idmv2 = 1
!             idv1 = 1
!             idv2 = 1
!             
!             if ((mstv1(1)-MSf(1,3,3))/MSf(1,3,3).gt.1.d-5) idmv1=0
!             if ((mstv2(1)-MSf(1,3,3))/MSf(1,3,3).gt.1.d-5) idmv2=0

!new version, added by Moritz to sort out potential NANs, which can appear e.g. inMQv1,MQv2 for special parameter-combinations (->negative argument in sqrt)

            idmv1 = 0
            idmv2 = 0
            idv1 = 1
            idv2 = 1
            
            if (dabs(mstv1(1)-MSf(1,3,3))/MSf(1,3,3).le.1.d-5) idmv1=1   !MO:should sort out potential NANs
            if (dabs(mstv2(1)-MSf(1,3,3))/MSf(1,3,3).le.1.d-5) idmv2=1
            
            stv1 = dabs(dreal(ustv1(1,1)*ustv1(1,2)))
            sbv1 = dabs(dreal(usbv1(1,1)*usbv1(1,2)))
            stv2 = dabs(dreal(ustv2(1,1)*ustv2(1,2)))
            sbv2 = dabs(dreal(usbv2(1,1)*usbv2(1,2)))

            sttv1 = dabs(dreal(ustv1(1,1)*ustv1(2,2) + ustv1(1,2)*ustv1(2,1)))
            sttv2 = dabs(dreal(ustv2(1,1)*ustv2(2,2) + ustv2(1,2)*ustv2(2,1)))

            if ((stv1.le.1d-4).or.(sttv1.le.1d-4).or.(sbv1.le.1d-4)) idv1 = 0
            if ((stv2.le.1d-4).or.(sttv2.le.1d-4).or.(sbv2.le.1d-4)) idv2 = 0

            iv1 = idv1*idmv1
            iv2 = idv2*idmv2

            if ((iv1.eq.0).and.(iv2.eq.0)) then
                    print*,'*** Warning: Both solutions incompatible ***'
                    renfail = 0

            else
*           solution 1 incompatible
                  if (iv1.eq.0) then
                        do m=1,2
                              MSf(m,3,3) = mstv2(m)
                              MSf2(m,3,3) = mstv22(m)
                              do n=1,2
                                USf(m,n,3,3) = ustv2(m,n)
                                USf(m,n,4,3) = usbv2(m,n)
                              enddo  
                        enddo      
*           solution 2 incompatible
                  else if (iv2.eq.0) then
                        do m=1,2
                              MSf(m,3,3) = mstv1(m)
                              MSf2(m,3,3) = mstv12(m)
                              do n=1,2
                                USf(m,n,3,3) = ustv1(m,n)
                                USf(m,n,4,3) = usbv1(m,n)
                              enddo  
                        enddo                                    
*           both solutions compatible
                  else
                     if ( dabs(mstv1(2)-MSf(2,3,3)).lt.dabs(mstv2(2)-MSf(2,3,3)) ) then
                           do m=1,2
                                 MSf(m,3,3) = mstv1(m)
                                 MSf2(m,3,3) = mstv12(m)
                                 do n=1,2
                                   USf(m,n,3,3) = ustv1(m,n)
                                   USf(m,n,4,3) = usbv1(m,n)
                                 enddo  
                           enddo                                       
                     else
                           do m=1,2
                                 MSf(m,3,3) = mstv2(m)
                                 MSf2(m,3,3) = mstv22(m)
                                 do n=1,2
                                   USf(m,n,3,3) = ustv2(m,n)
                                   USf(m,n,4,3) = usbv2(m,n)
                                 enddo  
                           enddo                                       
                     endif      
                  endif
            endif
            
            USf(2,1,3,3) = USf(2,1,3,3)*(-1)
            USf(2,2,3,3) = USf(2,2,3,3)*(-1)




! =====================================================
! ======= Mixed DR-OS Scheme No.2 - Counterterms ======
! =======                                        ======
! ============ mt,mb,Ab,At defined in DRbar ===========
! ========= mst1,msb1,msb2 defined on-shell ===========
! ========= theta, mst2 dependant parameters ==========
! ===================================================== 

*************************************************
* Alpha_s renormalization constant (DRbar, MSSM)
* deltaGs is in fact dGS/GS
         nf = 6
         DeltaGs = 0d0
         
! usual version         
         DeltaGs = (Alfas/(8d0*pi))*(UVdiv*(nf - 3*CV))
                  
! debugging version - splitting into SM and SUSY contrib

! ! SM contrib - gluon + ghost + quark loop
         DeltaGs = 0d0
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(2d0*nf/3d0 - 11d0*CV/3d0))
! ! SUSY contrib - gluino loop
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(1d0*nf/3d0 + 2d0*CV/3d0))

     
*************************************************
* Gluon wave-function renormalization constant

      call DMNLO_CalcdPiG(0d0)	        
      dZG = - dreal(dPiGT)     
   
! Ghost wave-function renormalization constant
      call DMNLO_CalcdPiC(0d0)
      dZC=-dreal(dPiC)
     
**************************            
* Quark mass counterterms

* DR counterterms for all quarks including the top quark
        do type=3,4 
          do gen=1,3
             dMf(type,gen) = Mf(type,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)
          enddo    
        enddo          
         
***********************************************
* Quark wave-function renormalization constants
	  do type=3,4
	    do gen=1,3
	        call DMNLO_CalcPif(Mf2(type,gen),type,gen)
	        call DMNLO_CalcdPif(Mf2(type,gen),type,gen)
	        dZfL(type,gen) = -PifVL(type,gen) - Mf2(type,gen)*(dPifVL(type,gen) + dPifVR(type,gen)) - Mf(type,gen)*dPifSL(type,gen) - Mf(type,gen)*dPifSR(type,gen)
	        dZfR(type,gen) = -PifVR(type,gen) - Mf2(type,gen)*(dPifVR(type,gen) + dPifVL(type,gen)) - Mf(type,gen)*dPifSR(type,gen) - Mf(type,gen)*dPifSL(type,gen)
	    enddo
	  enddo

***************************
* Squark mass counterterms

* OS for sbottom 1,2 and stop 1
	do type=3,4
	  do gen=1,3
	    do k=1,2
	
	      call DMNLO_CalcPiSf(MSf2(k,type,gen),type,gen)
        	dMSf2(k,type,gen) = PiSf(k,k,type,gen)
	
	    enddo
	  enddo
	enddo

* stop 2 dependant counterterm - defined BELOW

**************************************************
* Squark wave-function renormalization constants

	  do type=3,4
	    do gen=1,3
	      do k=1,2
	        do l=1,2
	          if (k.eq.l) then
	           call DMNLO_CalcdPiSf(MSf2(k,type,gen),type,gen)
	           dZSf(k, k, type,gen) = -dPiSf(k,k,type,gen)
	          else 
	           call DMNLO_CalcPiSf(MSf2(l,type,gen),type,gen)
	           dZSf(k, l, type,gen) = (2d0/(MSf2(k,type,gen) - MSf2(l,type,gen)))*PiSf(k,l,type,gen) 
	          endif
	        enddo
	      enddo
	    enddo
	  enddo

***************************
* Af counterterms in DRbar
      
         do type=3,4
           do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))

              dmfdr = Mf(type,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)            

              dthtm1m2(type,gen) = GS**2*CF/(16d0*pi**2)*
     &              ( (u11*u21-u22*u12)*((u11**2-u12**2)*MSf2(1,type,gen) + (u21**2-u22**2)*MSf2(2,type,gen)) + 
     &               4d0*(u11*u22+u12*u21)*MGl*Mf(type,gen))*UVdiv

              dmsf12 = GS**2*CF/(16d0*pi**2)*( (u11**2-u12**2)**2*MSf2(1,type,gen) + (u21*u11-u22*u12)**2*MSf2(2,type,gen) - 
     &                            MSf2(1,type,gen) + 8d0*(u11*u12)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv
    
              dmsf22 = GS**2*CF/(16d0*pi**2)*( (u21**2-u22**2)**2*MSf2(2,type,gen) + (u21*u11-u22*u12)**2*MSf2(1,type,gen) - 
     &                            MSf2(2,type,gen) + 8d0*(u21*u22)*Mf(type,gen)*MGl - 4d0*MGl**2 - 4d0*Mf(type,gen)**2)*UVdiv


              dM2LR = u11*u12*dmsf12 + u21*u22*dmsf22 + (u21*u12 + u11*u22)*dthtm1m2(type,gen)

              dAf(type,gen) = 1d0/Mf(type,gen)*(dM2LR - dmfdr/Mf(type,gen)*(u11*u12*MSf2(1,type,gen) + u21*u22*MSf2(2,type,gen)))

           enddo
         enddo
      
***************************************************
* Squark mixing angle & mixing matrix counterterms

         do type=3,4
           do gen=1,3
              u11=dreal(USf(1,1,type,gen))
              u12=dreal(USf(1,2,type,gen))
              u21=dreal(USf(2,1,type,gen))
              u22=dreal(USf(2,2,type,gen))
              
              if (type.eq.3) then                              
                  dM2LR = dMf(type,gen)*(Af(type,gen)-MUE/TB) + Mf(type,gen)*dAf(type,gen)
              else 
                  dM2LR = dMf(type,gen)*(Af(type,gen)-MUE*TB) + Mf(type,gen)*dAf(type,gen)
              endif
              
              dtht(type,gen) = 1d0/((u21*u12 + u11*u22)*(MSf2(1,type,gen)-MSf2(2,type,gen)))*(dM2LR - u11*u12*(dMSf2(1,type,gen) - dMSf2(2,type,gen)))
              dthtm1m2(type,gen) = 1d0/(u21*u12 + u11*u22)*(dM2LR - u11*u12*(dMSf2(1,type,gen) - dMSf2(2,type,gen)))
              
               do m=1,2
                do n=1,2
                   dUSf(m,n,4,3) = 0d0
                enddo
              enddo


              do m=1,2
                dUSf(1,1,type,gen) = dUSf(1,1,type,gen) + Epsij(1,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(1,2,type,gen) = dUSf(1,2,type,gen) + Epsij(1,m)*USf(m,2,type,gen)*dtht(type,gen)
                dUSf(2,1,type,gen) = dUSf(2,1,type,gen) + Epsij(2,m)*USf(m,1,type,gen)*dtht(type,gen)
                dUSf(2,2,type,gen) = dUSf(2,2,type,gen) + Epsij(2,m)*USf(m,2,type,gen)*dtht(type,gen)
              enddo
              
           enddo
         enddo

***************************************************
* Dependant counterterms dtheta, dUSf, dmst2

         ut11=dreal(USf(1,1,3,3))
         ut12=dreal(USf(1,2,3,3))
         ut21=dreal(USf(2,1,3,3))
         ut22=dreal(USf(2,2,3,3))
         
         ub11=dreal(USf(1,1,4,3))
         ub12=dreal(USf(1,2,4,3))
         ub21=dreal(USf(2,1,4,3))
         ub22=dreal(USf(2,2,4,3))
         
            
         dMSf2(2,3,3) = (ut21*ut12 + ut11*ut22)/(ut21*ut12)*
     &   (ub11**2*dMSf2(1,4,3) + ub21**2*dMSf2(2,4,3) + 2d0*ub11*ub21*dthtm1m2(4,3) - 2d0*Mf(4,3)*dMf(4,3)
     &   -ut11**2*dMSf2(1,3,3) + 2d0*Mf(3,3)*dMf(3,3)) -
     &   2d0*ut11/ut12*(dMf(3,3)*(Af(3,3)-MUE/TB) + Mf(3,3)*dAf(3,3) - ut11*ut12*dMSf2(1,3,3))          

         dtht(3,3) = 1d0/((ut21*ut12 + ut11*ut22)*(MSf2(1,3,3)-MSf2(2,3,3)))*
     &               (dMf(3,3)*(Af(3,3)-MUE/TB) + Mf(3,3)*dAf(3,3) - ut11*ut12*(dMSf2(1,3,3) - dMSf2(2,3,3)))

        do m=1,2
         do n=1,2
            dUSf(m,n,3,3) = 0d0
         enddo
       enddo

         do m=1,2
           dUSf(1,1,3,3) = dUSf(1,1,3,3) + Epsij(1,m)*USf(m,1,3,3)*dtht(3,3)
           dUSf(1,2,3,3) = dUSf(1,2,3,3) + Epsij(1,m)*USf(m,2,3,3)*dtht(3,3)
           dUSf(2,1,3,3) = dUSf(2,1,3,3) + Epsij(2,m)*USf(m,1,3,3)*dtht(3,3)
           dUSf(2,2,3,3) = dUSf(2,2,3,3) + Epsij(2,m)*USf(m,2,3,3)*dtht(3,3)
         enddo         
         

! ==========================================
!       only for debugging purposes
! ==========================================

       else if (renscheme.eq.99) then

! ==========================================
! Contributions for the gluons & quarks 
! ==========================================

*************************************************
* Alpha_s renormalization constant (DRbar, MSSM)
* deltaGs is in fact dGS/GS
         nf = 6
         DeltaGs = 0d0
         
! usual version         
         DeltaGs = (Alfas/(8d0*pi))*(UVdiv*(nf - 3*CV))
                  
! debugging version - splitting into SM and SUSY contrib

! ! SM contrib - gluon + ghost + quark loop
         DeltaGs = 0d0
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(2d0*nf/3d0 - 11d0*CV/3d0))
! ! SUSY contrib - gluino loop
         DeltaGs = DeltaGs + (Alfas/(8d0*pi))*(UVdiv*(1d0*nf/3d0 + 2d0*CV/3d0))



*************************************************
* Gluon wave-function renormalization constant

      call DMNLO_CalcdPiG(0d0)	        
      dZG = - dreal(dPiGT)
  
! Ghost wave-function renormalization constant
      call DMNLO_CalcdPiC(0d0)
      dZC=-dreal(dPiC)

	        	        
! ===== Mass renormalization constants =====
	do gen=1,3 
	  do type=3,4
	    call DMNLO_CalcPif(Mf2(type,gen),type,gen)
	    dMf(type,gen) = (Mf(type,gen)/2d0)*(PifVL(type,gen) + PifVR(type,gen)) + (PifSL(type,gen) + PifSR(type,gen))/2d0
	  enddo
	enddo

! ===== Wave-function renormalization constants =====
	do type=3,4
	  do gen=1,3
	        call DMNLO_CalcPif(Mf2(type,gen),type,gen)
	        call DMNLO_CalcdPif(Mf2(type,gen),type,gen)
 	        dZfL(type,gen) = -PifVL(type,gen) - Mf2(type,gen)*(dPifVL(type,gen) + dPifVR(type,gen)) - Mf(type,gen)*dPifSL(type,gen) - Mf(type,gen)*dPifSR(type,gen)
 	        dZfR(type,gen) = -PifVR(type,gen) - Mf2(type,gen)*(dPifVR(type,gen) + dPifVL(type,gen)) - Mf(type,gen)*dPifSR(type,gen) - Mf(type,gen)*dPifSL(type,gen)

	  enddo
	enddo

! ==========================================
! Contributions for the squarks 
! ==========================================

! ===== Mass renormalization constants =====

	do type=3,4
	  do gen=1,3
	    do k=1,2
		MSf2temp(k,type,gen) = MSf2(k,type,gen)	
	    enddo
	  enddo
	enddo
	
! ===== DRbar squark masses =====
!       (iteration)
	do i=1,3

	do type=3,4
	  do gen=1,3
	    do k=1,2

	      call DMNLO_CalcPiSf(MSf2temp(k,type,gen),type,gen)
	      MSf2temp(k,type,gen) = MSf2(k,type,gen) + PiSf(k,k,type,gen)

	    enddo
	  enddo
	enddo
	enddo


! ===== DRbar mass counterterm =====
	do type=3,4
	  do gen=1,3
	    do k=1,2
	          
	        call DMNLO_CalcPiSf(MSf2(k,type,gen),type,gen)
        	  dMSf2(k, type,gen) = PiSf(k,k,type,gen)
        	  
	    enddo
	  enddo
	enddo
	
! ===== Wave-function renormalization constants =====
	do type=3,4
	  do gen=1,3
	    do k=1,2
	      do l=1,2
	        if (k.eq.l) then
	         call DMNLO_CalcdPiSf(MSf2(k,type,gen),type,gen)
	         dZSf(k, k, type,gen) = -dPiSf(k,k,type,gen)
	        else 
	         call DMNLO_CalcPiSf(MSf2(l,type,gen),type,gen)
	         dZSf(k, l, type,gen) = (2d0/(MSf2(k,type,gen) - MSf2(l,type,gen)))*PiSf(k,l,type,gen) 
	        endif
	      enddo
	    enddo
	  enddo
	enddo

! ===== Mixing angle renormalization constants =====
	do type=3,4       
	  do gen=1,3
	    do k=1,2
	      do l=1,2
	        do m=1,2
	          dUSf(k, l, type,gen) = dUSf(k, l, type,gen) + 0.25D0*(dZSf(k, m, type,gen) - dconjg(dZSf(m, k, type,gen)))*USf(m, l, type,gen)
	        enddo
	      enddo
	    enddo
	  enddo
	enddo               


! Trilineare Kopplung

	do type=3,4
	  do gen=1,3

              dAf(type,gen) = (dMSf2(1,type,gen)*USf(1,1,type,gen)*USf(1,2,type,gen) +
     &                         dMSf2(2,type,gen)*USf(2,1,type,gen)*USf(2,2,type,gen) +
     &                         MSf(1,type,gen)**2*dUSf(1,1,type,gen)*USf(1,2,type,gen) +
     &                         MSf(1,type,gen)**2*USf(1,1,type,gen)*dUSf(1,2,type,gen) + 
     &                         MSf(2,type,gen)**2*dUSf(2,1,type,gen)*USf(2,2,type,gen) + 
     &                         MSf(2,type,gen)**2*USf(2,1,type,gen)*dUSf(2,2,type,gen) - 
     &  (MSf(1,type,gen)**2*USf(1,1,type,gen)*USf(1,2,type,gen) + 
     &   MSf(2,type,gen)**2*USf(2,1,type,gen)*USf(2,2,type,gen))*dMf(type,gen)/Mf(type,gen))/Mf(type,gen)
     
	  enddo
	enddo

        endif


999	return 
      end

! ===========================================      
! ===== QCD Self-Energies & derivatives =====
! ===========================================

! ============================================
! QCD derivative of Self-energy of the gluon 
! ============================================

	subroutine DMNLO_CalcdPiG(myP2)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	double precision myP2,ML0
	double complex g0,g1
	integer ttf,ttgen,k

	dPiGT = 0d0
        ML0 = 0d0
	
! ===== Squark loop =====	

	g0 = -GS
	g1 = -GS

	do ttf=3,4
		do ttgen=1,3
			do k=1,2
	
	    ML0 = MSf(k,ttf,ttgen)

	call myBget(myP2,ML0**2,ML0**2)
	call myBPget(myP2,ML0**2,ML0**2)

        dPiGT = dPiGT + TF*g0*g1*(3d0*xb0 - 12d0*ML0**2*xbp0)/(144d0*Pi**2)
        if (EpsPole.eq.0d0) dPiGT = dPiGT + TF*g0*g1*2d0/(144d0*Pi**2)
       
               enddo
            enddo
        enddo        

! ===== Quark loop =====	

	g0 = -GS
	g1 = -GS

	do ttf=3,4
           do ttgen=1,2
		
	        ML0 = 0d0


	call myBget(myP2,ML0**2,ML0**2)
	call myBPget(myP2,ML0**2,ML0**2)


        dPiGT = dPiGT + TF*g0*g1*(6d0*ML0**2*xbp0 + 3d0*xb0)/(36d0*Pi**2)
!        if (EpsPole.eq.0d0) dPiGT = dPiGT - TF*g0*g1/(36d0*Pi**2)

      
            enddo
        enddo
                
        do ttf=3,4
	    		
	ML0 = Mf(ttf,3)


	call myBget(myP2,ML0**2,ML0**2)
	call myBPget(myP2,ML0**2,ML0**2)


        dPiGT = dPiGT + TF*g0*g1*(6d0*ML0**2*xbp0 + 3d0*xb0)/(36d0*Pi**2)
!       if (EpsPole.eq.0d0) dPiGT = dPiGT - TF*g0*g1/(36d0*Pi**2)
        
      
        enddo        

        
! ===== Gluino loop =====	

	g0 = -GS
	g1 = -GS
		
	ML0 = MGl

	call myBget(myP2,ML0**2,ML0**2)
	call myBPget(myP2,ML0**2,ML0**2)

        dPiGT = dPiGT + CV*g0*g1*(6d0*ML0**2*xbp0 + 3d0*xb0)/(72d0*Pi**2)
        if (EpsPole.eq.0d0) dPiGT = dPiGT - CV*g0*g1/(72d0*Pi**2)
        

          
! ===== Gluon loop =====	

	g0 = -GS
	g1 = -GS
		
 	ML0 = 0d0
! 
 	call myBget(myP2,ML0**2,ML0**2)

        dPiGT = dPiGT + CV*g0*g1/(576d0*Pi**2)*(-57d0*xb0)
        if (EpsPole.eq.0d0) dPiGT = dPiGT + CV*g0*g1*67d0/(576d0*Pi**2)
        

! ===== Ghost loop =====
	
	g0 = -GS
	g1 = -GS
		
	ML0 = 0d0

	call myBget(myP2,ML0**2,ML0**2)

        dPiGT = dPiGT + CV*g0*g1/(576d0*Pi**2)*(-3d0*xb0)
        if (EpsPole.eq.0d0) dPiGT = dPiGT + CV*g0*g1/(576d0*Pi**2)
        
      
        end	

! =======================================
! QCD Self-energy of the quark  
! =======================================

	subroutine DMNLO_CalcPif(myP2,type,gen)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2

	PifVR(type,gen) = 0d0
	PifVL(type,gen) = 0d0
	PifSR(type,gen) = 0d0
	PifSL(type,gen) = 0d0

! ===== Gluon loop =====
	g1L = -GS
	g1R = -GS
	g2L = -GS
	g2R = -GS

	ML0 = 0d0 
	ML1 = Mf(type,gen)


	call myBget(myP2,ML0**2,ML1**2)

        PifVR(type,gen) = PifVR(type,gen) + CF/(16d0*pi**2)*2d0*g1R*g2R*(xb0 + xb1)
        PifVL(type,gen) = PifVL(type,gen) + CF/(16d0*pi**2)*2d0*g1L*g2L*(xb0 + xb1)
        PifSR(type,gen) = PifSR(type,gen) - CF/(16d0*pi**2)*4d0*g1R*g2L*ML1*(xb0)
        PifSL(type,gen) = PifSL(type,gen) - CF/(16d0*pi**2)*4d0*g1L*g2R*ML1*(xb0)



! ===== Gluino Loop =====

	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,type,gen)
	 	g1R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = -sqrt2*GS*USf(m,1,type,gen)

	 	ML0 = MSf(m,type,gen) 
	 	ML1 = MGl 

	 	call myBget(myP2,ML0**2,ML1**2)
         
		PifVR(type,gen) = PifVR(type,gen) + CF/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
         	PifVL(type,gen) = PifVL(type,gen) + CF/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
         	PifSR(type,gen) = PifSR(type,gen) + CF/(16d0*pi**2)*g1R*g2R*ML1*xb0
         	PifSL(type,gen) = PifSL(type,gen) + CF/(16d0*pi**2)*g1L*g2L*ML1*xb0

	enddo


      end	


! ===================================================
! QCD derivative of Self-energy of the quark  
! ===================================================

	subroutine DMNLO_CalcdPif(myP2,type,gen)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2


	dPifVR(type,gen) = 0d0
	dPifVL(type,gen) = 0d0
	dPifSR(type,gen) = 0d0
	dPifSL(type,gen) = 0d0

! ===== Gluon Loop =====
	g1L = -GS
	g1R = -GS
	g2L = -GS
	g2R = -GS

	ML0 = 0d0 
	ML1 = Mf(type,gen)

	call myBPget(myP2,ML1**2,ML0**2)
        dPifVR(type,gen) = dPifVR(type,gen) - CF/(16d0*pi**2)*2d0*g1R*g2R*xbp1
        dPifVL(type,gen) = dPifVL(type,gen) - CF/(16d0*pi**2)*2d0*g1L*g2L*xbp1
        dPifSR(type,gen) = dPifSR(type,gen) - CF/(16d0*pi**2)*4d0*g1R*g2L*ML1*xbp0
        dPifSL(type,gen) = dPifSL(type,gen) - CF/(16d0*pi**2)*4d0*g1L*g2R*ML1*xbp0

! ===== Gluino Loop =====

	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,type,gen)
	 	g1R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = -sqrt2*GS*USf(m,1,type,gen)

	 	ML0 = MSf(m,type,gen) 
	 	ML1 = MGl

	 	call myBPget(myP2,ML1**2,ML0**2)

         	dPifVR(type,gen) = dPifVR(type,gen) - CF/(16d0*pi**2)*g1R*g2L*xbp1
         	dPifVL(type,gen) = dPifVL(type,gen) - CF/(16d0*pi**2)*g1L*g2R*xbp1
         	dPifSR(type,gen) = dPifSR(type,gen) + CF/(16d0*pi**2)*g1R*g2R*ML1*xbp0
         	dPifSL(type,gen) = dPifSL(type,gen) + CF/(16d0*pi**2)*g1L*g2L*ML1*xbp0

	enddo

	end	
	
	
! =======================================
! QCD Self-energy of the gluino  
! =======================================

	subroutine DMNLO_CalcPigl(myP2)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l,nf
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2

	PiglVR = 0d0
	PiglVL = 0d0
	PiglSR = 0d0
	PiglSL = 0d0

!===== Gluon loop =====
	g1 = GS
	g2 = GS

	ML0 = 0d0 
	ML1 = MGl

	call myBget(myP2,ML0**2,ML1**2)

        PiglVR = PiglVR + Nc/(16d0*pi**2)*2d0*g1*g2*(xb0 + xb1)
        PiglVL = PiglVL + Nc/(16d0*pi**2)*2d0*g1*g2*(xb0 + xb1)
        
        PiglSR = PiglSR - Nc/(16d0*pi**2)*4d0*g1*g2*ML1*(xb0)
        PiglSL = PiglSL - Nc/(16d0*pi**2)*4d0*g1*g2*ML1*(xb0)

! ===== Squark/Quark Loop =====

	do m=1,2
	  do type=3,4
	    do gen=1,3

	 	g1R = -sqrt2*GS*USf(m,1,type,gen)
	 	g1L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = -sqrt2*GS*USf(m,1,type,gen)



	 	ML0 = MSf(m,type,gen) 
	 	ML1 = Mf(type,gen) 

	 	call myBget(myP2,ML0**2,ML1**2)
         
          PiglVR = PiglVR + 1/(32d0*pi**2)*g1R*g2L*(xb0 + xb1)
         	PiglVL = PiglVL + 1/(32d0*pi**2)*g1L*g2R*(xb0 + xb1)

          PiglSR = PiglSR + 1/(32d0*pi**2)*g1R*g2R*ML1*xb0
         	PiglSL = PiglSL + 1/(32d0*pi**2)*g1L*g2L*ML1*xb0


	    enddo
	  enddo
	enddo

! ===== Squark Antiquark Loop =====

	do m=1,2
	  do type=3,4
	    do gen=1,3

	 	g1L = -sqrt2*GS*USf(m,1,type,gen)
	 	g1R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = -sqrt2*GS*USf(m,1,type,gen)


	 	ML0 = MSf(m,type,gen) 
	 	ML1 = Mf(type,gen) 

	 	call myBget(myP2,ML0**2,ML1**2)
!          
          PiglVR = PiglVR - 1/(32d0*pi**2)*g2R*g1L*(xb0 + xb1)
         	PiglVL = PiglVL - 1/(32d0*pi**2)*g2L*g1R*(xb0 + xb1)

          PiglSR = PiglSR + 1/(32d0*pi**2)*g1R*g2R*ML1*xb0
         	PiglSL = PiglSL + 1/(32d0*pi**2)*g1L*g2L*ML1*xb0

	    enddo
	  enddo
	enddo

      end	


! ===================================================
! QCD derivative of Self-energy of the gluino
! ===================================================

	subroutine DMNLO_CalcdPigl(myP2)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2


	dPiglVR = 0d0
	dPiglVL = 0d0
	dPiglSR = 0d0
	dPiglSL = 0d0

! ===== Gluon Loop =====
	g1 = -GS
	g2 = -GS

	ML0 = 0d0 
	ML1 = MGl

	call myBPget(myP2,ML1**2,ML0**2)
	
	dPiglVR = dPiglVR - Nc/(16d0*pi**2)*2d0*g1*g2*xbp1
        dPiglVL = dPiglVL - Nc/(16d0*pi**2)*2d0*g1*g2*xbp1
        dPiglSR = dPiglSR - Nc/(16d0*pi**2)*4d0*g1*g2*ML1*xbp0
        dPiglSL = dPiglSL - Nc/(16d0*pi**2)*4d0*g1*g2*ML1*xbp0
        

! ===== Squark Quark Loop =====

	do m=1,2
	  do type=1,3
	    do gen=3,4

	 	g1L = -sqrt2*GS*USf(m,1,type,gen)
	 	g1R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = -sqrt2*GS*USf(m,1,type,gen)

	 	ML0 = MSf(m,type,gen) 
	 	ML1 = Mf(type,gen)

	 	call myBPget(myP2,ML1**2,ML0**2)

         	dPiglVR = dPiglVR - 1/(32d0*pi**2)*g1R*g2L*xbp1
         	dPiglVL = dPiglVL - 1/(32d0*pi**2)*g1L*g2R*xbp1
         	dPiglSR = dPiglSR - 1/(32d0*pi**2)*g1R*g2R*ML1*xbp0
         	dPiglSL = dPiglSL - 1/(32d0*pi**2)*g1L*g2L*ML1*xbp0


	    enddo
	  enddo
	enddo
	
	! ===== Squark Antiquark Loop =====

	do m=1,2
	  do type=1,3
	    do gen=3,4

	 	g1L = -sqrt2*GS*USf(m,1,type,gen)
	 	g1R = sqrt2*GS*USf(m,2,type,gen)
	 	g2L = sqrt2*GS*USf(m,2,type,gen)
	 	g2R = -sqrt2*GS*USf(m,1,type,gen)

	 	ML0 = MSf(m,type,gen) 
	 	ML1 = Mf(type,gen)

	 	call myBPget(myP2,ML1**2,ML0**2)

         	dPiglVL = dPiglVL - 1/(32d0*pi**2)*g1R*g2L*xbp1
         	dPiglVR = dPiglVR - 1/(32d0*pi**2)*g1L*g2R*xbp1
         	dPiglSL = dPiglSL - 1/(32d0*pi**2)*g1R*g2R*ML1*xbp0
         	dPiglSR = dPiglSR - 1/(32d0*pi**2)*g1L*g2L*ML1*xbp0


	    enddo
	  enddo
	enddo

	end		
! =========================================
! QCD Self-energy of the squark 
! =========================================

	subroutine DMNLO_CalcPiSf(myP2,type,gen)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2

	complex*16 B02,B02p
	double precision A01, Krondelta

	do k=1,2
	  do l=1,2
	    PiSf(k,l,type,gen) = 0d0
 	  enddo
	enddo

	do k=1,2
	  do l=1,2

! ===== Gluon Loop =====
	    g1 = -GS*Krondelta(k,l)
	    g2 = -GS*Krondelta(k,l)

	    ML0 = 0d0 
	    ML1 = MSf(k,type,gen)

            PiSf(k,l,type,gen) = PiSf(k,l,type,gen) - CF/(16d0*pi**2)*g1*g2*(-A01(ML1) + 2d0*(ML1**2 + myP2)*B02(myP2,ML0,ML1))

! ===== Gluino Loop =====
	    g1L = sqrt2*GS*USf(l,2,type,gen)
	    g1R = -sqrt2*GS*USf(l,1,type,gen)
	    g2L = -sqrt2*GS*USf(k,1,type,gen)
	    g2R = sqrt2*GS*USf(k,2,type,gen)

	    ML0 = Mf(type,gen) 
	    ML1 = MGl

            PiSf(k,l,type,gen) = PiSf(k,l,type,gen) - CF/(16d0*pi**2)*(2d0*(g1L*g2L + g1R*g2R)*ML0*ML1*B02(myP2,ML0,ML1) + 
     &       	(g1R*g2L + g1L*g2R)*(A01(ML0) + A01(ML1) + (ML0**2 + ML1**2 - myP2)*B02(myP2,ML0,ML1)))

! ===== Squark Loop =====
	    do m=1,2

	   	g1 = -GS**2*(USf(m,1,type,gen)*USf(l,1,type,gen)-USf(m,2,type,gen)*USf(l,2,type,gen))*
     &                (USf(m,1,type,gen)*USf(k,1,type,gen)-USf(m,2,type,gen)*USf(k,2,type,gen))

	   	ML0 = MSf(m,type,gen)

           	PiSf(k,l,type,gen) = PiSf(k,l,type,gen) - CF/(16d0*pi**2)*g1*A01(ML0)

	    enddo

	  enddo
	enddo

	end

! ====================================================
! QCD derivative of Self-energy of the squark 
! ====================================================

	subroutine DMNLO_CalcdPiSf(myP2,type,gen)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer type,gen,k,m,l
	double precision myP2
	double complex g1L,g1R,g2L,g2R,ML0,ML1,g1,g2

	complex*16 B02,B02p
	double precision A01, Krondelta

	do k=1,2
	  do l=1,2
	    dPiSf(k,l,type,gen) = 0d0 
 	  enddo
	enddo

	do k=1,2
	  do l=1,2

! ===== Gluon Loop =====
	    g1 = -GS*Krondelta(k,l)
	    g2 = -GS*Krondelta(k,l)

	    ML0 = 0d0 
	    ML1 = MSf(k,type,gen)

            dPiSf(k,l,type,gen) = dPiSf(k,l,type,gen) - CF/(16d0*pi**2)*g1*g2*(2d0*(ML1**2 + myP2)*B02p(myP2,ML0,ML1) + 2d0*B02(myP2,ML0,ML1))

! ===== Gluino Loop =====
	    g1L = sqrt2*GS*USf(l,2,type,gen)
	    g1R = -sqrt2*GS*USf(l,1,type,gen)
	    g2L = -sqrt2*GS*USf(k,1,type,gen)
	    g2R = sqrt2*GS*USf(k,2,type,gen)

	    ML0 = Mf(type,gen) 
	    ML1 = MGl

            dPiSf(k,l,type,gen) = dPiSf(k,l,type,gen) - CF/(16d0*pi**2)*(2d0*(g1L*g2L + g1R*g2R)*ML0*ML1*B02p(myP2,ML0,ML1) + 
     &       	(g1R*g2L + g1L*g2R)*((ML0**2 + ML1**2 - myP2)*B02p(myP2,ML0,ML1) - B02(myP2,ML0,ML1)))

	  enddo
	enddo

	end


! ==================================================================
! QCD Self-energy of the bottom quark for threshold & top iteration
! ==================================================================

	subroutine DMNLO_CalcPifth(myP2,Q,t,g)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer t,g,k,m,l
	double precision myP2,Q,ML0,ML1
	double precision UVtemp,Qtemp
	double complex g1L,g1R,g2L,g2R,g1,g2

* store global divergence & scale
      UVtemp = UVdiv
      Qtemp = muSc
      UVdiv = 0d0
      muSc = Q

	PifVRth(t,g) = 0d0
	PifVLth(t,g) = 0d0
	PifSRth(t,g) = 0d0
	PifSLth(t,g) = 0d0

      if ((t.eq.3).and.(g.eq.3)) then
! ===== Gluon loop =====
      	g1L = -GS
      	g1R = -GS
      	g2L = -GS
      	g2R = -GS
      
      	ML0 = 0d0 
      	ML1 = dsqrt(myP2)
            
      	call myBget(myP2,ML0**2,ML1**2)
      
              PifVRth(t,g) = PifVRth(t,g) + CF/(16d0*pi**2)*2d0*g1R*g2R*(xb0 + xb1)
              PifVLth(t,g) = PifVLth(t,g) + CF/(16d0*pi**2)*2d0*g1L*g2L*(xb0 + xb1)
              PifSRth(t,g) = PifSRth(t,g) - CF/(16d0*pi**2)*4d0*g1R*g2L*ML1*(xb0)
              PifSLth(t,g) = PifSLth(t,g) - CF/(16d0*pi**2)*4d0*g1L*g2R*ML1*(xb0)

      end if

! ===== Gluino-Squark Loop =====

	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,t,g)
	 	g1R = sqrt2*GS*USf(m,2,t,g)
	 	g2L = sqrt2*GS*USf(m,2,t,g)
	 	g2R = -sqrt2*GS*USf(m,1,t,g)

	 	ML0 = MSf(m,t,g) 
	 	ML1 = MGl 

	 	call myBget(myP2,ML0**2,ML1**2)

         	PifVRth(t,g) = PifVRth(t,g) + CF/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
        	PifVLth(t,g) = PifVLth(t,g) + CF/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
         	PifSRth(t,g) = PifSRth(t,g) + CF/(16d0*pi**2)*g1R*g2R*ML1*xb0
         	PifSLth(t,g) = PifSLth(t,g) + CF/(16d0*pi**2)*g1L*g2L*ML1*xb0

	enddo

      if ((t.eq.4).and.(g.eq.3)) then
! ===== Chargino-Stop Loop =====

	do m=1,2
	    do k=1,2

	 	g1L = -EL/SW*VCha(k,1)*USf(m,1,3,3) + EL*Mf(3,3)/(sqrt2*MW*SW*SB)*VCha(k,2)*USf(m,2,3,3) 
	 	g1R = EL*dsqrt(myP2)/(sqrt2*MW*SW*CB)*UCha(k,2)*USf(m,1,3,3) 
	 	g2L = EL*dsqrt(myP2)/(sqrt2*MW*SW*CB)*UCha(k,2)*USf(m,1,3,3) 
	 	g2R = -EL/SW*VCha(k,1)*USf(m,1,3,3) + EL*Mf(3,3)/(sqrt2*MW*SW*SB)*VCha(k,2)*USf(m,2,3,3) 

	 	ML0 = MSf(m,3,3) 
	 	ML1 = MCha(k) 

	 	call myBget(myP2,ML0**2,ML1**2)

         	PifVRth(t,g) = PifVRth(t,g) + 1d0/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
         	PifVLth(t,g) = PifVLth(t,g) + 1d0/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
         	PifSRth(t,g) = PifSRth(t,g) + 1d0/(16d0*pi**2)*g1R*g2R*ML1*xb0
         	PifSLth(t,g) = PifSLth(t,g) + 1d0/(16d0*pi**2)*g1L*g2L*ML1*xb0

	    enddo
	enddo

      end if

* revert global divergence to original value
      UVdiv = UVtemp
      muSc = Qtemp

	end	

! ==================================================================
! QCD Self-energy of the bottom quark for testing
! ==================================================================

	subroutine DMNLO_CalcPiftest(myP2,Q,t,g,iloop)

#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

	integer t,g,k,m,l,iloop
	double precision myP2,Q,ML0,ML1
	double precision UVtemp,Qtemp
	double complex g1L,g1R,g2L,g2R,g1,g2

* store global divergence & scale
      UVtemp = UVdiv
      Qtemp = muSc
      UVdiv = 0d0
      muSc = Q

	PifVRtest(t,g) = 0d0
	PifVLtest(t,g) = 0d0
	PifSRtest(t,g) = 0d0
	PifSLtest(t,g) = 0d0

      if ((iloop.eq.0).or.(iloop.eq.1)) then
! ===== Gluon loop =====
      	g1L = -GS
      	g1R = -GS
      	g2L = -GS
      	g2R = -GS
      
      	ML0 = 0d0 
      	ML1 = dsqrt(myP2)
            
      	call myBget(myP2,ML0**2,ML1**2)
      
              PifVRtest(t,g) = PifVRtest(t,g) + CF/(16d0*pi**2)*2d0*g1R*g2R*(xb0 + xb1)
              PifVLtest(t,g) = PifVLtest(t,g) + CF/(16d0*pi**2)*2d0*g1L*g2L*(xb0 + xb1)
              PifSRtest(t,g) = PifSRtest(t,g) - CF/(16d0*pi**2)*4d0*g1R*g2L*ML1*(xb0)
              PifSLtest(t,g) = PifSLtest(t,g) - CF/(16d0*pi**2)*4d0*g1L*g2R*ML1*(xb0)

      end if

! ===== Gluino-Squark Loop =====

      if ((iloop.eq.0).or.(iloop.eq.2)) then
	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,t,g)
	 	g1R = sqrt2*GS*USf(m,2,t,g)
	 	g2L = sqrt2*GS*USf(m,2,t,g)
	 	g2R = -sqrt2*GS*USf(m,1,t,g)

	 	ML0 = MSf(m,t,g) 
	 	ML1 = MGl 

	 	call myBget(myP2,ML0**2,ML1**2)

         	PifVRtest(t,g) = PifVRtest(t,g) + CF/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
        	PifVLtest(t,g) = PifVLtest(t,g) + CF/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
         	PifSRtest(t,g) = PifSRtest(t,g) + CF/(16d0*pi**2)*g1R*g2R*ML1*xb0
         	PifSLtest(t,g) = PifSLtest(t,g) + CF/(16d0*pi**2)*g1L*g2L*ML1*xb0

	enddo

      if ((t.eq.4).and.(g.eq.3)) then
! ===== Chargino-Stop Loop =====

	do m=1,2
	    do k=1,2

	 	g1L = -EL/SW*VCha(k,1)*USf(m,1,3,3) + EL*Mf(3,3)/(sqrt2*MW*SW*SB)*VCha(k,2)*USf(m,2,3,3) 
	 	g1R = EL*dsqrt(myP2)/(sqrt2*MW*SW*CB)*UCha(k,2)*USf(m,1,3,3) 
	 	g2L = EL*dsqrt(myP2)/(sqrt2*MW*SW*CB)*UCha(k,2)*USf(m,1,3,3) 
	 	g2R = -EL/SW*VCha(k,1)*USf(m,1,3,3) + EL*Mf(3,3)/(sqrt2*MW*SW*SB)*VCha(k,2)*USf(m,2,3,3) 

	 	ML0 = MSf(m,3,3) 
	 	ML1 = MCha(k) 

	 	call myBget(myP2,ML0**2,ML1**2)

         	PifVRtest(t,g) = PifVRtest(t,g) + 1d0/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
         	PifVLtest(t,g) = PifVLtest(t,g) + 1d0/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
         	PifSRtest(t,g) = PifSRtest(t,g) + 1d0/(16d0*pi**2)*g1R*g2R*ML1*xb0
         	PifSLtest(t,g) = PifSLtest(t,g) + 1d0/(16d0*pi**2)*g1L*g2L*ML1*xb0

	    enddo
	enddo

      end if
      end if
* revert global divergence to original value
      UVdiv = UVtemp
      muSc = Qtemp
	
      end	      
  
! ===================================================
! QCD derivative of Self-energy of ghosts 
! ===================================================

      subroutine DMNLO_CalcdPiC(myP2)
      implicit none
#include "../DMNLO_Model.h"
#include "../DMNLO_RenScheme.h"
#include "../DMNLO_LoopIntegrals.h"

      double precision myP2
  
      call myBget(myP2,0d0,0d0)
  
      dPiC=-GS**2*CV/(32*Pi**2)*xb0
  
      if (EpsPole.eq.0d0) dPiC = dPiC+GS**2*CV/(32*Pi**2)
  
      end 
  