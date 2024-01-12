! ===============================================================================================
! Subroutine to set the kinematics of the ingoing 4-momenta and magnitude of in/out-going momenta

	subroutine QQ2xx_SetKinematics(mytree)

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

      integer mytree,a,swap,charge1,charge2,charge3,charge4,charge5

! lepton-masses to regularize 2->3-processes with leptonic final state

        ME = 0.000510998928d0 
        ME2 = ME**2  
        Mf(2,1) = ME   
        Mf2(2,1) = ME2       
        MM = 0.1056583715d0 
        MM2 = MM**2
        Mf(2,2) = MM
        Mf2(2,2) = MM2

!on-shell lepton masses used in 2->3 for first two generations (remaining fermion masses already defined in DMNLO_MSSMinit.F)

	Mfos(1,1) = 0
	Mfos(2,1) = ME
	Mfos(1,2) = 0
	Mfos(2,2) = MM

!Vector-Masses

        MV(1) = 0d0             !Gluon
        MV(2) = 0d0             !Photon
        MV(3) = MZ              !Z-Boson
        MV(4) = MW              !W_Boson


! ========== Reminder : indices convention =======================================================	   
c sftype1 :    Initial Sfermion1 type (1 : Sneutrino, 2 : Slepton, 3 : up Squark, 4 : down Squark) 
c sfgen1 :     Initial Sfermion1 generation (1, 2, 3)
c sfmeigen1 :  Initial Sfermion1 mass eigenstate (1, 2 : L, R or 1, 2)
c sftype2 :    Initial Sfermion2 type (1 : Sneutrino, 2 : Slepton, 3 : up Squark, 4 : down Squark) 
c sfgen2 :     Initial Sfermion2 generation (1, 2, 3)
c sfmeigen2 :  Initial Sfermion2 mass eigenstate (1, 2 : L, R or 1, 2)
c f3type :     Final fermion 3 type (1 : neutrino, 2 : lepton, 3 : up quark, 4 : down quark) 
c f3gen :      Final fermion 3 generation (1, 2, 3)
c f4type :     Final fermion 4 type (1 : neutrino, 2 : lepton, 3 : up quark, 4 : down quark) 
c f4gen :      Final fermion 4 generation (1, 2, 3)
c xtype :      Final boson 3 type (1 : gluon, 2 : photon, 3 : Z, 4 : W+, 5: h, 6: H, 7: A, 8: H+)
c xhtype :     Final boson 3 type (1 : h, 2 : H, 3 : A0, 4 : G0, 5 : H+, 6 : G+)                      !needed for some couplings
c ytype :      Final boson 4 type (1 : gluon, 2 : photon, 3 : Z, 4 : W-, 5: h, 6: H, 7: A, 8: H-)
c yhtype :     Final boson 4 type (1 : h, 2 : H, 3 : A0, 4 : G0, 5 : H-, 6 : G-)                      !needed for some couplings

!Convention: 1 is a sQuark (momentum p1), 2 is an anti-squark (momentum p2), 3 is a boson(+ charge) or a fermion (momentum k1), 4 is a boson(- charge) or an anti-fermion (momentum k2), 5 is the Gluon in 2->3-processes (momentum k3)

c ========== Example : stop1 sTop1 -> W+ H- process ==========
c	sftype = 3
c	sfgen = 3
c	sfmeigen = 1
c       xtype = 4
c       ytype = 8
c       f3type = 0
c       f3gen = 0
c       f4type = 0
c       f4gen = 0
c       charge1 = 2
c       charge2 = -2
c       charge3 = 3
c       charge4 = -3
c ========== Example : stop1 sTop1 -> t T process ============
c	sftype = 3
c	sfgen = 3
c	sfmeigen = 1
c       f3type = 3
c       f3gen = 3
c       f4type = 4
c       f4gen = 4
c       xtype = 0
c       ytype = 0
c       charge1 = 2
c       charge2 = -2
c       charge3 = 2
c       charge4 = -2
! ================================================================================================

! ===== Index setting ============================
! here indices are set in a convenient form for the use with the micrOMEGAs-improve-cross-section-function and the QQ2xx-convention:
! 1 is sparticle, 2 is anti-sparticle, 3 is boson(+) or fermion, 4 is boson(-) or antifermion 
! if n1 and n2 or n3 and n4 are both negative or positive, these special-cases will be picked up by tree-fct-call etc.

        if(n1<0) then

             swap=n1
             n1=n2
             n2=swap

        endif

        if(n3<0) then

             swap=n3
             n3=n4
             n4=swap

        endif

        if((n4.eq.24).or.(n4.eq.37)) then

             swap=n3
             n3=n4
             n4=swap

        endif
        
! ========== initial state particle #1 - sfermions ========== 

	select case (abs(n1))            

! --- Squarks ---                
        case (1000001)
                sftype1 =  4
                sfgen1 = 1
                sfmeigen1 = 1
                charge1 = -1
        case (1000002)
                sftype1 =  3
                sfgen1 = 1
                sfmeigen1 = 1
                charge1 = 2
        case (1000003)
                sftype1 =  4
                sfgen1 = 2
                sfmeigen1 = 1
                charge1 = -1
        case (1000004)
                sftype1 =  3
                sfgen1 = 2
                sfmeigen1 = 1
                charge1 = 2
        case (1000005)
                sftype1 =  4
                sfgen1 = 3
                sfmeigen1 = 1
                charge1 = -1
        case (1000006)
                sftype1 =  3
                sfgen1 = 3
                sfmeigen1 = 1
                charge1 = 2
        case (2000001)
                sftype1 =  4
                sfgen1 = 1
                sfmeigen1 = 2
                charge1 = -1
        case (2000002)
                sftype1 =  3
                sfgen1 = 1
                sfmeigen1 = 2
                charge1 = 2
        case (2000003)
                sftype1 =  4
                sfgen1 = 2
                sfmeigen1 = 2
                charge1 = -1
        case (2000004)
                sftype1 =  3
                sfgen1 = 2
                sfmeigen1 = 2
                charge1 = 2
        case (2000005)
                sftype1 =  4
                sfgen1 = 3
                sfmeigen1 = 2
                charge1 = -1
        case (2000006)
                sftype1 =  3
                sfgen1 = 3
                sfmeigen1 = 2
                charge1 = 2
                
! --- Sleptons ---
        case (1000011)
                sftype1 =  2
                sfgen1 = 1
                sfmeigen1 = 1
                charge1 = -3
        case (1000012)
                sftype1 =  1
                sfgen1 = 1
                sfmeigen1 = 1
                charge1 = 0
        case (1000013)
                sftype1 =  2
                sfgen1 = 2
                sfmeigen1 = 1
                charge1 = -3
        case (1000014)
                sftype1 =  1
                sfgen1 = 2
                sfmeigen1 = 1
                charge1 = 0
        case (1000015)
                sftype1 =  2
                sfgen1 = 3
                sfmeigen1 = 1
                charge1 = -3
        case (1000016)
                sftype1 =  1
                sfgen1 = 3
                sfmeigen1 = 1
                charge1 = 0
        case (2000011)
                sftype1 =  2
                sfgen1 = 1
                sfmeigen1 = 2
                charge1 = -3
        case (2000013)
                sftype1 =  2
                sfgen1 = 2
                sfmeigen1 = 2
                charge1 = -3
        case (2000015)
                sftype1 =  2
                sfgen1 = 3
                sfmeigen1 = 2
                charge1 = -3
	case DEFAULT
		sftype1 = 0
		sfgen1 = 0
		sfmeigen1 = 0 
                charge1 = 111
        
        end select
   
! ========== initial state particle #2 - sfermions ========== 

	select case (abs(n2))

! --- Squarks ---
        case (1000001)
                sftype2 =  4
                sfgen2 = 1
                sfmeigen2 = 1
                charge2 = 1
        case (1000002)
                sftype2 =  3
                sfgen2 = 1
                sfmeigen2 = 1
                charge2 = -2
        case (1000003)
                sftype2 =  4
                sfgen2 = 2
                sfmeigen2 = 1
                charge2 = 1
        case (1000004)
                sftype2 =  3
                sfgen2 = 2
                sfmeigen2 = 1
                charge2 = -2
        case (1000005)
                sftype2 =  4
                sfgen2 = 3
                sfmeigen2 = 1
                charge2 = 1
        case (1000006)
                sftype2 =  3
                sfgen2 = 3
                sfmeigen2 = 1  
                charge2 = -2              
        case (2000001)
                sftype2 =  4
                sfgen2 = 1
                sfmeigen2 = 2
                charge2 = 1
        case (2000002)
                sftype2 =  3
                sfgen2 = 1
                sfmeigen2 = 2
                charge2 = -2
        case (2000003)
                sftype2 =  4
                sfgen2 = 2
                sfmeigen2 = 2
                charge2 = 1
        case (2000004)
                sftype2 =  3
                sfgen2 = 2
                sfmeigen2 = 2
                charge2 = -2
        case (2000005)
                sftype2 =  4
                sfgen2 = 3
                sfmeigen2 = 2
                charge2 = 1
        case (2000006)
                sftype2 =  3
                sfgen2 = 3
                sfmeigen2 = 2
                charge2 = -2
               
! --- Sleptons ---
        case (1000011)
                sftype2 =  2
                sfgen2 = 1
                sfmeigen2 = 1
                charge2 = 3
        case (1000012)
                sftype2 =  1
                sfgen2 = 1
                sfmeigen2 = 1
                charge2 = 0
        case (1000013)
                sftype2 =  2
                sfgen2 = 2
                sfmeigen2 = 1
                charge2 = 3
        case (1000014)
                sftype2 =  1
                sfgen2 = 2
                sfmeigen2 = 1
                charge2 = 0
        case (1000015)
                sftype2 =  2
                sfgen2 = 3
                sfmeigen2 = 1
                charge2 = 3
        case (1000016)
                sftype2 =  1
                sfgen2 = 3
                sfmeigen2 = 1
                charge2 = 0
        case (2000011)
                sftype2 =  2
                sfgen2 = 1
                sfmeigen2 = 2
                charge2 = 3
        case (2000013)
                sftype2 =  2
                sfgen2 = 2
                sfmeigen2 = 2
                charge2 = 3
        case (2000015)
                sftype2 =  2
                sfgen2 = 3
                sfmeigen2 = 2
                charge2 = 3
	case DEFAULT
		sftype2 = 0
		sfgen2 = 0
		sfmeigen2 = 0   
                charge2 = 222             
                
        end select
	
! ========== final state particle #3 - boson ========== 	

	select case (abs(n3))

! gluon final state	
	case (21)                                                !since I negelect flavor-mixing I set propgen=sfgen1 etc.
		xtype = 1
                charge3 = 0
		
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
		
! photon final state
	case (22)
		xtype = 2
                charge3 = 0	
	
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
		
! Z final state
	case (23)
		xtype = 3
                charge3 = 0
	
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
		
! W+ final state
	case (24)
		xtype = 4
                charge3 = 3
		
                if ((sftype1.eq.1).or.(sftype1.eq.3)) then       !Type-Index for propagator masses
                proptypeu=sftype1
                proptypeuc=sftype1
                proptypet=sftype1+1
                proptypetc=sftype1+1
                else
                proptypeu=sftype1
                proptypeuc=sftype1
                proptypet=sftype1-1
                proptypetc=sftype1-1
                endif

                propgenu=sfgen1
                propgenuc=sfgen1
                propgent=sfgen1
                propgentc=sfgen1
		
! h final state
	case (25)
		xtype = 5	
                xhtype = 1
                charge3 = 0
	
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
		
! H final state
	case (35)
		xtype = 6
                xhtype = 2
                charge3 = 0
	
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
			
! H3 final state
	case (36)
		xtype = 7
                xhtype = 3
                charge3 = 0
	
                proptypeu=sftype1
                proptypeuc=sftype1
                propgenu=sfgen1
                propgenuc=sfgen1
                proptypet=sftype1
                proptypetc=sftype1
                propgent=sfgen1
                propgentc=sfgen1
				
! H+ final state
	case (37)
		xtype = 8
                xhtype = 5 
                charge3 = 3

                if ((sftype1.eq.1).or.(sftype1.eq.3)) then       !Type-Index for propagatormasses
                proptypeu=sftype1
                proptypeuc=sftype1
                proptypet=sftype1+1
                proptypetc=sftype1+1
                else
                proptypeu=sftype1
                proptypeuc=sftype1
                proptypet=sftype1-1
                proptypetc=sftype1-1
                endif

                propgenu=sfgen1
                propgenuc=sfgen1
                propgent=sfgen1
                propgentc=sfgen1

	case DEFAULT
		xtype = 0
                xhtype = 0
                proptypeu=0
                proptypeuc=0
                propgenu=0
                propgenuc=0
                proptypet=0
                proptypetc=0
                propgent=0
                propgentc=0
                charge3 = 333

	end select

! ========== final state particle #3 - fermion ========== 	

	select case (abs(n3))

! --- Quarks ---
        case (1)
                f3type =  4
                f3gen = 1
                charge3 = -1
        case (2)
                f3type =  3
                f3gen = 1
                charge3 = 2
        case (3)
                f3type =  4
                f3gen = 2
                charge3 = -1
        case (4)
                f3type =  3
                f3gen = 2
                charge3 = 2
        case (5)
                f3type =  4
                f3gen = 3
                charge3 = -1
        case (6)
                f3type =  3
                f3gen = 3  
                charge3 = 2              

! --- Leptons ---
        case (11)
                f3type =  2
                f3gen = 1
                charge3 = -3
        case (12)
                f3type =  1
                f3gen = 1
                charge3 = 0
        case (13)
                f3type =  2
                f3gen = 2
                charge3 = -3
        case (14)
                f3type =  1
                f3gen = 2
                charge3 = 0
        case (15)
                f3type =  2
                f3gen = 3
                charge3 = -3
        case (16)
                f3type =  1
                f3gen = 3
                charge3 = 0
	case DEFAULT
		f3type = 0
          	f3gen = 0
                      
        end select        

! ========== final state particle #4 - boson ========== 

	select case (abs(n4))
	
! gluon final state	
	case (21)
		ytype = 1
                charge4 = 0
				
! photon final state
	case (22)
		ytype = 2
                charge4 = 0
				
! Z final state
	case (23)
		ytype = 3
                charge4 = 0
				
! W- final state
	case (24)
		ytype = 4
                charge4 = -3
		
                if ((sftype1.eq.1).or.(sftype1.eq.3)) then       !Type-Index for propagatormasses
                proptypeu=sftype1+1
                proptypeuc=sftype1+1
                else
                proptypeu=sftype1-1
                proptypeuc=sftype1-1
                endif

                propgenu=sfgen1
                propgenuc=sfgen1
	
! h final state
	case (25)
		ytype = 5		
                yhtype = 1
                charge4 = 0                
		
! H final state
	case (35)
		ytype = 6
                yhtype = 2
                charge4 = 0
			
! H3 final state
	case (36)
		ytype = 7
                yhtype = 3
                charge4 = 0
			
! H- final state
	case (37)
		ytype = 8
                yhtype = 5
                charge4 = -3

                if ((sftype1.eq.1).or.(sftype1.eq.3)) then       !Type-Index for propagatormasses
                proptypeu=sftype1+1
                proptypeuc=sftype1+1
                else
                proptypeu=sftype1-1
                proptypeuc=sftype1-1
                endif

                propgenu=sfgen1
                propgenuc=sfgen1
		
	case DEFAULT
		ytype = 0
                yhtype = 0
                proptypeu=0
                proptypeuc=0
                propgenu=0
                propgenuc=0
                proptypet=0
                proptypetc=0
                propgent=0
                propgentc=0
                charge4 = 444
		
	end select

! ========== final state particle #4 - fermion ========== 	

	select case (abs(n4))

! --- Quarks ---
        case (1)
                f4type =  4
                f4gen = 1
                charge4 = 1
        case (2)
                f4type =  3
                f4gen = 1
                charge4 = -2
        case (3)
                f4type =  4
                f4gen = 2
                charge4 = 1
        case (4)
                f4type =  3
                f4gen = 2
                charge4 = -2
        case (5)
                f4type =  4
                f4gen = 3
                charge4 = 1
        case (6)
                f4type =  3
                f4gen = 3    
                charge4 = -2            

! --- Leptons ---
        case (11)
                f4type =  2
                f4gen = 1
                charge4 = 3
        case (12)
                f4type =  1
                f4gen = 1
                charge4 = 0
        case (13)
                f4type =  2
                f4gen = 2
                charge4 = 3
        case (14)
                f4type =  1
                f4gen = 2
                charge4 = 0
        case (15)
                f4type =  2
                f4gen = 3
                charge4 = 3
        case (16)
                f4type =  1
                f4gen = 3
                charge4 = 0
	case DEFAULT
		f4type = 0
          	f4gen = 0
                      
        end select

! ========== final state particle #5 - boson ========== 

	select case (abs(n5))
	
! gluon final state	
	case (21)

		ztype = 1
                charge5 = 0
				
! photon final state
	case (22)
		ztype = 2
                charge5 = 0
		
! Z final state
	case (23)
		ztype = 3
                charge5 = 0	
! W+- final state
	case (24)
		ztype = 4
                charge5 = 3
    		
	case DEFAULT
		ztype = 0
                charge5 = 555
               		
	end select

! === Set charges right (particle <-> antiparticle)===
!if the incoming or outgoing particles are particle-antiparticle-pairs, this part doesn't change anything,
!in the remaining cases the charges still need to be set correctly-> done here

        if(n1<0) charge1=-charge1
        if(n2>0) charge2=-charge2
        if(n3<0) charge3=-charge3
        if(n4>0) charge4=-charge4
        if(n5<0) charge5=-charge5
        
! === total charges (for checks on charge conservation) ===

        total_chargeout = charge3 + charge4              !needed for some couplings/amplitudes (VHgen/HVgen)
        total_charge22 = charge1+charge2-total_chargeout !check on charge conservation
        total_charge23 = total_charge22-charge5          !check on charge conservation

! === Summation-limits for (s,t,u)-channel-exchange-particles ===

        if(sftype1.eq.sftype2) then                      
  
        sumstartH = 1
        sumendH = 4
        sumstartV = 2
        sumstartVV = 1
        sumendV = 3

        else

        sumstartH = 5
        sumendH = 6
        sumstartV = 4
        sumstartVV = 4
        sumendV = 4

        endif

! === Check processes and set mass parameters ===

	if ((sftype1.eq.0).or.(sftype2.eq.0).or.((f3type.eq.0).and.(xtype.eq.0)).or.((f4type.eq.0).and.(ytype.eq.0)).or.((f3type.eq.0).and.(ytype.eq.0)).or.((f4type.eq.0).and.(xtype.eq.0)).or.(sftype1.lt.3).or.(sftype2.lt.3)
     &      .or.((sfgen1.ne.sfgen2).and.(f3type.eq.0).and.(f4type.eq.0))) then         !case of incoming sLeptons not covered,no generation-mixing in initial and final state, the last part shall set 
                                                                                       !processeses with flavormixing including Higgs- or Vector final states to zero (ugly but necessary, because VVSfSf/HHSfSf would otherwise 
                                                                                       !give nonzero contributions although not written for the flavor mixing case)
                                                                                                                                                                                                          
	  m1 = 0d0
	  m2 = 0d0
	  m3 = 0d0
	  m4 = 0d0
          m5 = 0d0 

!           print*,"  "
!           print *,"-------- this case is not covered by QQ2xx --------"
!           print*,"  "

          abort22 = 3
          abort23 = 1

        else
           
!2->2-part

          if(total_charge22.ne.0) then          
          
!             print*,"  "
!             print*,'total_charge22 not 0', total_charge22
!             print*,"  "

            abort22 = 3

          else if((mytree.eq.0).and.((xtype.eq.1).or.(ytype.eq.1).or.(f3type.gt.2).or.(f4type.gt.2))) then

!             print*,"  "
!             print *,"-------- strongly interacting final-states exist only at tree-level --------"
!             print*,"  "

            if((xtype.eq.1).and.(ytype.gt.1)) then           !Gluon-elw. Boson -> needed for/in on-shell subtraction

              abort22 = 2

            else if((xtype.gt.1).and.(ytype.eq.1)) then      !elw. Boson-Gluon -> needed for/in on-shell subtraction

              abort22 = 2

            else
              
              abort22 = 1

            endif

          else

            abort22 = 0

          endif

!2->3-part

          if(abort22.ne.0) then

            abort23 = 1

          else if(total_charge23.ne.0) then

!             print*,"  "
!             print*,'total_charge23 not 0', total_charge23
!             print*,"  "

            abort23 = 1

          else if((mytree.eq.0).and.(ztype.eq.0)) then

!             print*,"  "
!             print *,"-------- this 2->3-choice is not covered by QQ2xx --------"
!             print*,"  "

            abort23 = 1

          elseif((mytree.eq.0).and.((ztype.eq.2).or.(ztype.eq.3).or.(ztype.eq.4))) then

!             print*,"  "
! 	      print *,"-------- this 2->3-choice is only partly covered by QQ2xx --------"
!             print*,"  "

            abort23 = 1

          else

            abort23 = 0

          endif


          m1 = MSf(sfmeigen1, sftype1, sfgen1)                  
	
          m2 = MSf(sfmeigen2, sftype2, sfgen2)

	
          if((xtype.eq.0).and.(ytype.eq.0)) then                   

	    m3 = Mf(f3type,f3gen)                                  
	
            m4 = Mf(f4type,f4gen)
         
          else                                                     

            set_third_mass: select case (xtype)

	    case(1) set_third_mass 
	   	m3 = 0 
  	    case(2) set_third_mass
  	   	m3 = 0 
	    case(3) set_third_mass
	   	m3 = MZ
	    case(4) set_third_mass 
	   	m3 = MW
	    case(5) set_third_mass
	   	m3 = Mh0
	    case(6) set_third_mass
	   	m3 = MHH
	    case(7) set_third_mass
	   	m3 = MA0
	    case(8) set_third_mass
	   	m3 = MHp
	    case DEFAULT
		m3 = 0d0
	    end select set_third_mass

	    set_fourth_mass: select case (ytype)

	    case(1) set_fourth_mass 
	    	m4 = 0 
  	    case(2) set_fourth_mass
  	   	m4 = 0 
	    case(3) set_fourth_mass
	   	m4 = MZ
	    case(4) set_fourth_mass 
	   	m4 = MW
	    case(5) set_fourth_mass
	   	m4 = Mh0
	    case(6) set_fourth_mass
	   	m4 = MHH
	    case(7) set_fourth_mass
	   	m4 = MA0
	    case(8) set_fourth_mass
	   	m4 = MHp

	    case DEFAULT
		m4 = 0d0
	    end select set_fourth_mass

          endif

          set_fifth_mass: select case (ztype)

	  case(1) set_fifth_mass 
	   	m5 = 0 
  	  case(2) set_fifth_mass
  	   	m5 = 0 
	  case(3) set_fifth_mass
	   	m5 = MZ
          case(4) set_fifth_mass
	   	m5 = MW
	  case DEFAULT
		m5 = 0d0
	  end select set_fifth_mass
           
	endif
	
	end
