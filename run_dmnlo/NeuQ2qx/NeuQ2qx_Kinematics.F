! ===============================================================================================
! Subroutine to set the kinematics of the ingoing 4-momenta and magnitude of in/out-going momenta

		subroutine NeuQ2qx_SetKinematics()
!		subroutine NeuQ2qx_SetKinematics(n1,n2,n3,n4)

    	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      integer ntemp



! ========== Reminder : indices convention =======================================================	
c chimeigen : Gaugino index (1,2,3,4 : Neutralino 1,2,3,4 and 5,6 : Chargino 1,2)       
c sftype :    Initial Sfermion type (1 : Sneutrino, 2 : Slepton, 3 : up Squark, 4 : down Squark) 
c sfgen :     Initial Sfermion generation (1, 2, 3)
c sfmeigen :  Initial Sfermion mass eigenstate (1, 2 : L, R or 1, 2)
c ftype :     Final fermion type (1 : neutrino, 2 : lepton, 3 : up quark, 4 : down quark) 
c fgen :      Final fermion generation (1, 2, 3)
c xtype :     Final boson type (1 : gluon, 2 : photon, 3 : Z, 4 : W, 5: h, 6: H, 7: A, 8: H+)
c ========== Example : top gluon process ==========
c	chimeigen = 1 
c	sftype = 3
c	sfgen = 3
c	sfmeigen = 1
c	ftype = 3
c	fgen = 3
c	xtype = 1
! ================================================================================================



! ===== Index setting ============================
! here indices are set
! because of later used "selectProcess" order is guaranteed

! sort particles so that n1 corresponds to neutralinos / charginos and n2 to sfermions 
	if(((n2.ge.1000022).and.(n2.le.1000025)).or.((n2.eq.1000035).and.(n2.eq.1000037))) then
		ntemp = n2
		n2 = n1
		n1 = ntemp
	endif
! sort particles so that n3 corresponds to leptons and n4 to SM bosons
	if(((n4.ge.1).and.(n4.le.6)).or.((n4.ge.11).and.(n4.le.16))) then
		ntemp = n3
		n3 = n4
		n4 = ntemp
	endif
! ========== initial state particle #1 - neutralinos / charginos ==========        

	select case (n1)
	case (1000022)
	        chimeigen = 1
        case (1000023)
	        chimeigen = 2
        case (1000025)
	        chimeigen = 3
        case (1000035)
	        chimeigen = 4
        case (1000024)
	        chimeigen = 5
        case (1000037)
	        chimeigen = 6
	case DEFAULT
		chimeigen = 0
        end select
	
        
! ========== initial state particle #2 - sfermions ========== 

	select case (n2)

! --- Squarks ---

        case (1000001)
                sftype =  4
                sfgen = 1
                sfmeigen = 1
        case (1000002)
                sftype =  3
                sfgen = 1
                sfmeigen = 1
        case (1000003)
                sftype =  4
                sfgen = 2
                sfmeigen = 1
        case (1000004)
                sftype =  3
                sfgen = 2
                sfmeigen = 1
        case (1000005)
                sftype =  4
                sfgen = 3
                sfmeigen = 1
        case (1000006)
                sftype =  3
                sfgen = 3
                sfmeigen = 1
                
        case (2000001)
                sftype =  4
                sfgen = 1
                sfmeigen = 2
        case (2000002)
                sftype =  3
                sfgen = 1
                sfmeigen = 2
        case (2000003)
                sftype =  4
                sfgen = 2
                sfmeigen = 2
        case (2000004)
                sftype =  3
                sfgen = 2
                sfmeigen = 2
        case (2000005)
                sftype =  4
                sfgen = 3
                sfmeigen = 2
        case (2000006)
                sftype =  3
                sfgen = 3
                sfmeigen = 2
                
! --- Sleptons ---

        case (1000011)
                sftype =  2
                sfgen = 1
                sfmeigen = 1
        case (1000012)
                sftype =  1
                sfgen = 1
                sfmeigen = 1
        case (1000013)
                sftype =  2
                sfgen = 2
                sfmeigen = 1
        case (1000014)
                sftype =  1
                sfgen = 2
                sfmeigen = 1
        case (1000015)
                sftype =  2
                sfgen = 3
                sfmeigen = 1
        case (1000016)
                sftype =  1
                sfgen = 3
                sfmeigen = 1
        case (2000011)
                sftype =  2
                sfgen = 1
                sfmeigen = 2
        case (2000013)
                sftype =  2
                sfgen = 2
                sfmeigen = 2
        case (2000015)
                sftype =  2
                sfgen = 3
                sfmeigen = 2

	case DEFAULT
		sftype = 0
		sfgen = 0
		sfmeigen = 0
                
        end select
	
! ========== final state particle #3 - fermions ========== 	

	select case (n3)

! --- Quarks ---

        case (1)
                ftype =  4
                fgen = 1
        case (2)
                ftype =  3
                fgen = 1
        case (3)
                ftype =  4
                fgen = 2
        case (4)
                ftype =  3
                fgen = 2
        case (5)
                ftype =  4
                fgen = 3
        case (6)
                ftype =  3
                fgen = 3                

! --- Sleptons ---

        case (11)
                ftype =  2
                fgen = 1
        case (12)
                ftype =  1
                fgen = 1
        case (13)
                ftype =  2
                fgen = 2
        case (14)
                ftype =  1
                fgen = 2
        case (15)
                ftype =  2
                fgen = 3
        case (16)
                ftype =  1
                fgen = 3
	case DEFAULT
		ftype = 0
		fgen = 0
                      
        end select
        
        
! ========== final state particle #4 - bosons ========== 
! at the moment only the neutralino as incoming particle considered!
! in t- / u- channel the first particle of all possible exchange particles is as default value set

	select case (n4)
	
! gluon final state	
	case (21)
		xtype = 1
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
! incommented line - why fixed on stop1 in t-channel? somehow necessary for general tree level?		
		ptmeigen = sfmeigen
!		ptmeigen = 1

		pumeigen = 1
		
! photon final state
	case (22)
		xtype = 2
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
c		ptmeigen = sfmeigen
		ptmeigen = 1

		pumeigen = 1
		
! Z final state
	case (23)
		xtype = 3
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
		ptmeigen = 1
		
		pumeigen = 1
		
! W+ final state
	case (24)
		xtype = 4
		
		pstype = sftype
		psgen = sfgen
		
		if ((sftype.eq.3).or.(sftype.eq.4)) then
		pttype = mod(sftype-2,2)+3
		else 
		pttype = mod(sftype,2)+1
		endif
		ptgen = sfgen
		ptmeigen = 1
		
		pumeigen = 5		
		
! h final state
	case (25)
		xtype = 5
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
		ptmeigen = 1
		
		pumeigen = 1
		
! H final state
	case (35)
		xtype = 6
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
		ptmeigen = 1
		
		pumeigen = 1
		
! H3 final state
	case (36)
		xtype = 7
		
		pstype = sftype
		psgen = sfgen
		
		pttype = sftype
		ptgen = sfgen
		ptmeigen = mod(sfmeigen,2)+1
		
		pumeigen = 1
		
! H+ final state
	case (37)
		xtype = 8
		
		pstype = sftype
		psgen = sfgen
		
		if ((sftype.eq.3).or.(sftype.eq.4)) then
		pttype = mod(sftype-2,2)+3
		else 
		pttype = mod(sftype,2)+1
		endif
		ptgen = sfgen
		ptmeigen = 1
		
		pumeigen = 5

	case DEFAULT
		xtype = 0
		
	end select

	
! ===== Set mass parameters =====

	if ((chimeigen.eq.0).or.(sftype.eq.0).or.(ftype.eq.0).or.(xtype.eq.0)) then 

	    m1 = 0d0
	    m2 = 0d0
	    m3 = 0d0
	    m4 = 0d0
	    print *,"-------- this case is not covered by NeuQ2qx --------"
		stop
	else

	if (chimeigen.le.5d0) then
		m1 = MNeu(chimeigen)
	else
		m1 = MCha(chimeigen)
	endif
	
	m2 = MSf(sfmeigen, sftype, sfgen)
	
	m3 = Mf(ftype, fgen)
	
	set_fourth_mass: select case (xtype)
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
	
	end