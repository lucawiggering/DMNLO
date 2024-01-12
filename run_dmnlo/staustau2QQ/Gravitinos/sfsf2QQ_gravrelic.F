 !**************************************************************************************************************************
	  program relicofgravitinos
	 
!!  calculates the gravitino relic density, dependent on its mass and the TR 
!!  temperature after inflation. It uses the stau relic density for a particular scenario
!!  where staus are NLSPs. 



     !REFERENCES:
     ![1] PRADLER: Thermal Gravitino Production and Collider Tests of Leptogenesis (arXiv: 0608344, 2007)
     ![2] LINDE/STAROBINSKY: Reheating after inflation (1995)
     ![3] STRUMIA/RYCHKOV: Thermal production of gravitinos (arXiv: 0701104, 2007)
     ![4] POSPELOV: Constraints on SUSY models from catalytic primordial nucleosynthesis (arXiv: 0807.4287, 2008)
     
        implicit none

      integer x,y
      double precision grelicMO,grelicLO,grelicNLO
      double precision TR,staurelicMO
      double precision staurelicLO,staurelicNLO,staumass,gluinomass
      double precision ratio,ierror
      double precision gmass
      
      !scheme1 
      
!       staurelicMO = 0.447
!       staurelicLO = 0.385
!       staurelicNLO = 0.402
!       
      !scheme2  (we checked dependence on renormalisation scheme)
      
      staurelicMO = 0.447
      staurelicLO = 0.427
      staurelicNLO = 0.402
      
      
      
      staumass = 1808.    !can be modified for other NLSPs
      gluinomass = 5000

      gmass =10.   !gravitino mass, undetermined yet
      TR = 10.**7

      open(unit=20,file='sfsf2QQ_gravBJ_NLOscheme2.txt',status='new')
 
      do x=0,500
       do y=0,750
    
      gmass = 300.+0.4*y   !!mass of gravitino with phenomenol. constraints
      TR = 10.**7+0.1*x*10.**6  !!TR temperatur in the early, post-inflationary universe, see [2]
      
      !constraint: formulae below only valid for T_R > 10⁷ 
      
      grelicMO = (gluinomass/1000d0)**2*TR/(10.**10)*(100./gmass)     !thermal production, see [3]
     -    +    staurelicMO*gmass/staumass                             ! non-thermal production, see [1]
     
      grelicLO = (gluinomass/1000d0)**2*TR/(10.**10)*(100./gmass)
     -    +    staurelicLO*gmass/staumass
     
      grelicNLO = (gluinomass/1000d0)**2*TR/(10.**10)*(100./gmass)
     -    +    staurelicNLO*gmass/staumass
     
      ratio = 100*(staurelicNLO*gmass/staumass)/                     !a dominant non-th. production is important, otherwise our NLO corrections 
     -    ((gluinomass/1000d0)**2*TR/(10.**10)*(100./gmass)+         !do not affect the relic density any more  
     -    staurelicNLO*gmass/staumass)
     
           
      write(20,*) gmass,TR,grelicMO,grelicLO,grelicNLO,ratio 
      
      !output: gravitino-LSP relic density, dependent on their mass and the reheating temperature T_R
      
      ! be aware of lifetime constraints for NLSP, see e.g. [4]
      
	enddo       
	enddo 
	
	close(unit=20)  
 

      end program
