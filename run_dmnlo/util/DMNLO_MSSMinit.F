! DMNLO_MSSMinit.F

      subroutine DMNLO_ModelPara(MOswitch,filename)
      implicit none

#include "DMNLO_Model.h"
#include "SLHA.h"
#include "PDG.h"

      integer error
      double complex slhadata(nslhadata)
      double precision Halpha,MB2OS,MBDR,mMS2mMS,alsmb,AS 
      character*(*) filename ! name of the SLHA file
      integer MOswitch


! ===== READ SLHA FILE =====
      call SLHARead(error, slhadata, filename, 1)
      if( error .ne. 0 ) stop "SLHA read error"
		 
! ===== Gaugino mass parameters =====
      M_1 = MSoft_M1
      M_2 = MSoft_M2
      MGl = Mass_MGl

! ===== Higgs sector parameters =====
      MUE = HMix_MUE
      TB  = HMix_TB

! ===== Slepton mass parameters =====
      MSL(1) = MSoft_MSL(1) 
      MSL(2) = MSoft_MSL(2) 
      MSL(3) = MSoft_MSL(3) 
      MSE(1) = MSoft_MSE(1)
      MSE(2) = MSoft_MSE(2)
      MSE(3) = MSoft_MSE(3)
	  
! ===== Trilinear couplings =====
      Af(2,1) = Af_Af(1,1,2)
      Af(2,2) = Af_Af(2,2,2)
      Af(2,3) = Af_Af(3,3,2)

! ===== Squark mass parameters =====
      MSQ(1) = MSoft_MSQ(1)
      MSQ(2) = MSoft_MSQ(2)
      MSQ(3) = MSoft_MSQ(3)
      MSU(1) = MSoft_MSU(1) 
      MSU(2) = MSoft_MSU(2) 
      MSU(3) = MSoft_MSU(3) 
      MSD(1) = MSoft_MSD(1)
      MSD(2) = MSoft_MSD(2)
      MSD(3) = MSoft_MSD(3)

! ===== Trilinear couplings =====
      Af(3,1) = Af_Af(1,1,3)
      Af(4,1) = Af_Af(1,1,4)
      Af(3,2) = Af_Af(2,2,3)
      Af(4,2) = Af_Af(2,2,4)
      Af(3,3) = Af_Af(3,3,3)
      Af(4,3) = Af_Af(3,3,4)

      Xf(2,1) = Af(2,1) - dconjg(MUE)*TB
      Xf(3,1) = Af(3,1) - dconjg(MUE)/TB
      Xf(4,1) = Af(4,1) - dconjg(MUE)*TB

      Xf(2,2) = Af(2,2) - dconjg(MUE)*TB
      Xf(3,2) = Af(3,2) - dconjg(MUE)/TB
      Xf(4,2) = Af(4,2) - dconjg(MUE)*TB

      Xf(2,3) = Af(2,3) - dconjg(MUE)*TB
      Xf(3,3) = Af(3,3) - dconjg(MUE)/TB
      Xf(4,3) = Af(4,3) - dconjg(MUE)*TB

! ===== Gauge boson and ghost masses (for Feynman gauge) =====
      MZ = SMInputs_MZ
      MZ2 = MZ**2
      wZ = 2.4952d0  ! PDG 2022 value

      MG0 = MZ
      MG02 = MG0**2
      wG0 = wZ

      MW = Mass_MW
      MW2 = MW**2
      wW = 2.085d0 ! PDG 2022 value

      MGp = MW
      MGp2 = MGp**2
      wGp = wW

! ===== Electroweak mixing =====	  
      CW = MW/MZ
      CW2 = CW**2
      SW2 = 1-CW2
      SW = dsqrt(SW2)

      if(MOswitch.eq.1) then
           SW = 0.481d0
           SW2 = SW**2
           CW2 = 1d0 - SW**2
           CW = dsqrt(CW2)
           MW = CW*MZ
      endif

      C2W = CW2-SW2            
      S2W = 2*SW*CW            
      TW = SW/CW
      TW2 = TW**2

! ===== Electroweak coupling =====	  
      Alfa = 1d0/SMInputs_invAlfaMZ
      Alfa2 = Alfa**2

! ===== Strong coupling =====	  
      AlfasMZ = SMInputs_AlfasMZ

! ===== Lepton masses =====
      ME = 0d0
      MM = 0d0
      ML = SMInputs_Mtau 

! ===== Quark masses (continued) =====
      MU = 0d0
      MD = 0d0
      MS = 0d0
      MC = 0d0

! MTos is mt_pole
	MTos = SMInputs_Mt
      MT = MTos

      ME2=ME**2
      MM2=MM**2
      ML2=ML**2

      MD2=MD**2
      MS2=MS**2
      MB2=MB**2

      MU2=MU**2
      MC2=MC**2
      MT2=MT**2

! ===== Higgs masses and widths =====
      Mh0 = Mass_Mh0    
      MHH = Mass_MHH
      MA0 = Mass_MA0  
      MHp = Mass_MHp 

      wh0 = SLHAGetDecay(slhadata, PDG_h0, 0)    ! h0
      wHH = SLHAGetDecay(slhadata, PDG_HH, 0)    ! H0
      wA0 = SLHAGetDecay(slhadata, PDG_A0, 0)    ! A0 
      wHp = SLHAGetDecay(slhadata, PDG_Hp, 0)    ! H+ 

      Mh02 = Mh0**2
      MHH2 = MHH**2
      MA02 = MA0**2
      MHp2 = MHp**2

! ===== Effective Higgs mixing parameters =====
      Halpha = Alpha_Alpha
      SA = dsin(Halpha)
      CA = dcos(Halpha)
	
! ===== Different bottom masses =====
! MBSM is mb^MS_SM (mb)
      MBSM = SMInputs_Mb
! MB is mb^MS_SM (MZ)
      alsmb =  4d0*Pi*AS(MBSM, MZ**2, AlfasMZ/(4d0*Pi), 5, 3, 5)  ! alpha_s in the MSbar-scheme with 5 active quark flavors at the scale mb(mb) 
      MB = mMS2mMS(MBSM,alsmb,AlfasMZ, 5, 4)
! MBos is the pole mass
      MBos = MB2OS(MBSM,MBSM,0)


      Mfm(1) = MB
      Mfm(2) = MB
      Mfm(3) = MB
      Mfm(5) = MB
      Mfm(4) = Mfm(5)

      end

! ======================================================================
      subroutine DMNLO_ModelIni(filename)
      implicit none
#include "DMNLO_Model.h"
#include "SLHA.h"
#include "PDG.h"

      integer i, j, t, g,error 
      double complex slhadata(nslhadata)
      character*(*) filename



      call SLHARead(error, slhadata, filename, 1)
      if( error .ne. 0 ) stop "SLHA read error"


! ===== CKM matrix (diagonal) =====
      CKM(1,1) = 1
      CKM(2,2) = 1
      CKM(3,3) = 1
      CKM(2,1) = 0
      CKM(3,1) = 0
      CKM(1,2) = 0
      CKM(3,2) = 0
      CKM(1,3) = 0
      CKM(2,3) = 0

! ===== Fermion masses =====
      Mf(1,1) = 0
      Mf(2,1) = ME
      Mf(3,1) = MU
      Mf(4,1) = MD
      Mf(1,2) = 0
      Mf(2,2) = MM
      Mf(3,2) = MC
      Mf(4,2) = MS
      Mf(1,3) = 0
      Mf(2,3) = ML
      Mf(3,3) = MT
      Mf(4,3) = MB

      Mf2(1,1) = 0
      Mf2(2,1) = ME2
      Mf2(3,1) = MU2
      Mf2(4,1) = MD2
      Mf2(1,2) = 0
      Mf2(2,2) = MM2
      Mf2(3,2) = MC2
      Mf2(4,2) = MS2
      Mf2(1,3) = 0
      Mf2(2,3) = ML2
      Mf2(3,3) = MT2
      Mf2(4,3) = MB2
	
! ===== Fermion masses correlated with Higgs =====
! h
      MfH(1,1,5) = 0
      MfH(2,1,5) = ME
      MfH(3,1,5) = MU
      MfH(4,1,5) = MB
      MfH(1,2,5) = 0
      MfH(2,2,5) = MM
      MfH(3,2,5) = MC
      MfH(4,2,5) = MS
      MfH(1,3,5) = 0
      MfH(2,3,5) = ML
      MfH(3,3,5) = MTscale
      MfH(4,3,5) = MB
! H	
      MfH(1,1,6) = 0
      MfH(2,1,6) = ME
      MfH(3,1,6) = MU
      MfH(4,1,6) = 0d0
      MfH(1,2,6) = 0
      MfH(2,2,6) = MM
      MfH(3,2,6) = MC
      MfH(4,2,6) = 0d0
      MfH(1,3,6) = 0d0
      MfH(2,3,6) = ML
      MfH(3,3,6) = MTscale
      MfH(4,3,6) = MB
! H3	
      MfH(1,1,7) = 0
      MfH(2,1,7) = ME
      MfH(3,1,7) = MU
      MfH(4,1,7) = 0d0
      MfH(1,2,7) = 0
      MfH(2,2,7) = MM
      MfH(3,2,7) = MC
      MfH(4,2,7) = 0d0
      MfH(1,3,7) = 0
      MfH(2,3,7) = ML
      MfH(3,3,7) = MTscale
      MfH(4,3,7) = MB
	
! ===== on-shell fermion masses =====	
      Mfos(1,1) = 0
      Mfos(2,1) = ME
      Mfos(3,1) = MU
      Mfos(4,1) = MD
      Mfos(1,2) = 0
      Mfos(2,2) = MM
      Mfos(3,2) = MC
      Mfos(4,2) = 0
      Mfos(1,3) = 0
      Mfos(2,3) = ML
      Mfos(3,3) = MTos
      Mfos(4,3) = MBos

! ===== Fermion widths =====
	
      width(1) = 0
      width(2) = 0
      width(3) = SLHAGetDecay(slhadata, PDG_top, 0)
      width(4) = 0

! ===== sFermion widths =====

      WSf(1,1,1) = 0d0                          !elw. sparticle-widths so far not needed
      WSf(1,2,1) = 0d0
      WSf(2,2,1) = 0d0

      WSf(1,1,2) = 0d0
      WSf(1,2,2) = 0d0
      WSf(2,2,2) = 0d0

      WSf(1,1,3) = 0d0
      WSf(1,2,3) = 0d0
      WSf(2,2,3) = 0d0

      WSf(1,3,1) = SLHAGetDecay(slhadata, PDG_sup1, 0)  ! ~uL
      WSf(1,4,1) = SLHAGetDecay(slhadata, PDG_sdown1, 0)  ! ~dL
      WSf(2,3,1) = SLHAGetDecay(slhadata, PDG_sup2, 0)  ! ~uR
      WSf(2,4,1) = SLHAGetDecay(slhadata, PDG_sdown2, 0)  ! ~dR
      WSf(1,3,2) = SLHAGetDecay(slhadata, PDG_scharm1, 0)  ! ~cL
      WSf(1,4,2) = SLHAGetDecay(slhadata, PDG_sstrange1, 0)  ! ~sL
      WSf(2,3,2) = SLHAGetDecay(slhadata, PDG_scharm2, 0)  ! ~cR
      WSf(2,4,2) = SLHAGetDecay(slhadata, PDG_sstrange2, 0)  ! ~sR
      WSf(1,3,3) = SLHAGetDecay(slhadata, PDG_stop1, 0)  ! ~t1

      WSf(1,4,3) = SLHAGetDecay(slhadata, PDG_sbottom1, 0)  ! ~b1
      WSf(2,3,3) = SLHAGetDecay(slhadata, PDG_stop2, 0)  ! ~t2
      WSf(2,4,3) = SLHAGetDecay(slhadata, PDG_sbottom2, 0)  ! ~b2

      do i=1,2                                !sometimes MO-width equals zero (e.g. for stop1) if possible decays kinematically highly surpressed->NANs in 2->3-processes due to unregularized propagators
      do j=3,4                              !nicht schön,Ergebniss zeigt auch zum Teil eine gewisse Abh. von Breite (d.h. für Breite 0.01->0,001 variiert WQ-Ergebnis teilweise um 20%)
      do t= 1,3

      if(WSf(i,j,t).lt.0.01) then
            WSf(i,j,t) = 0.01
      endif

      enddo
      enddo
      enddo

	wst(1) = 0d0
	wst(2) = 0d0

! ===== Electroweak parameters =====
	EL = sqrt(4*pi*Alfa)
	SW = sqrt(SW2)
	
	g2h1 = EL/SW
	g2h2 = g2h1**2
	

! ===== Higgs mixing =====
	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = 1 - CB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = CB2 - SB2

	CA2 = CA**2
	SA2 = SA**2
	S2A = 2*SA*CA
	C2A = CA2 - SA2
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA

! ===== Coupling constants =====
	Alfas2MZ = AlfasMZ**2
	GSMZ = sqrt(4*pi*AlfasMZ)

! ===== Gluino mass =====
	MGl2 = MGl**2

! ===== Calculate masses here if desired =====
      if(imass.eq.1) then
            call DMNLO_SfermionMasses
            call DMNLO_CharginoMasses
            call DMNLO_NeutralinoMasses
      else

! ===== Sfermion masses and mixing =====

      MSf(1,1,1) = Mass_MSf(1,1,1)
      MSf(1,2,1) = Mass_MSf(1,2,1)
      MSf(2,2,1) = Mass_MSf(2,2,1)

      MSf(1,1,2) = Mass_MSf(1,1,2)
      MSf(1,2,2) = Mass_MSf(1,2,2)
      MSf(2,2,2) = Mass_MSf(2,2,2)

      MSf(1,1,3) = Mass_MSf(1,1,3)
      MSf(1,2,3) = Mass_MSf(1,2,3)
      MSf(2,2,3) = Mass_MSf(2,2,3)

      MSf(1,3,1) = Mass_MSf(1,3,1)
      MSf(1,4,1) = Mass_MSf(1,4,1)
      MSf(2,3,1) = Mass_MSf(2,3,1)
      MSf(2,4,1) = Mass_MSf(2,4,1)

      MSf(1,3,2) = Mass_MSf(1,3,2)
      MSf(1,4,2) = Mass_MSf(1,4,2)
      MSf(2,3,2) = Mass_MSf(2,3,2)
      MSf(2,4,2) = Mass_MSf(2,4,2)

      MSf(1,3,3) = Mass_MSf(1,3,3)
      MSf(1,4,3) = Mass_MSf(1,4,3)
      MSf(2,3,3) = Mass_MSf(2,3,3)
      MSf(2,4,3) = Mass_MSf(2,4,3)

      do i=1,2
      do t=1,4
      do g=1,3

      MSf2(i,t,g) = MSf(i,t,g)**2

      enddo
      enddo
      enddo

      do t=1,4
      do g=1,3

      USf(1,1,t,g) = 1d0
      USf(1,2,t,g) = 0d0
      USf(2,1,t,g) = 0d0
      USf(2,2,t,g) = 1d0

      enddo
      enddo
 
	
      USf(1,1,2,3) = SfMix_USf(1,1,2)
      USf(1,2,2,3) = SfMix_USf(1,2,2)
      USf(2,1,2,3) = SfMix_USf(2,1,2)
      USf(2,2,2,3) = SfMix_USf(2,2,2)

      USf(1,1,3,3) = SfMix_USf(1,1,3) 
      USf(1,2,3,3) = SfMix_USf(1,2,3) 
      USf(2,1,3,3) = SfMix_USf(2,1,3) 
      USf(2,2,3,3) = SfMix_USf(2,2,3) 

      USf(1,1,4,3) = SfMix_USf(1,1,4) 
      USf(1,2,4,3) = SfMix_USf(1,2,4) 
      USf(2,1,4,3) = SfMix_USf(2,1,4) 
      USf(2,2,4,3) = SfMix_USf(2,2,4) 
	
! ===== Neutralino masses and mixing =====

      MNeu(1) = Mass_MNeu(1) 
      MNeu(2) = Mass_MNeu(2) 
      MNeu(3) = Mass_MNeu(3) 
      MNeu(4) = Mass_MNeu(4) 

      do i=1,4
            MNeu2(i) = MNeu(i)**2 
      enddo

      ZNeu(1,1) = NMix_ZNeu(1,1)
      ZNeu(1,2) = NMix_ZNeu(1,2)
      ZNeu(1,3) = NMix_ZNeu(1,3)
      ZNeu(1,4) = NMix_ZNeu(1,4)
      ZNeu(2,1) = NMix_ZNeu(2,1)
      ZNeu(2,2) = NMix_ZNeu(2,2)
      ZNeu(2,3) = NMix_ZNeu(2,3)
      ZNeu(2,4) = NMix_ZNeu(2,4)
      ZNeu(3,1) = NMix_ZNeu(3,1)
      ZNeu(3,2) = NMix_ZNeu(3,2)
      ZNeu(3,3) = NMix_ZNeu(3,3)
      ZNeu(3,4) = NMix_ZNeu(3,4)
      ZNeu(4,1) = NMix_ZNeu(4,1)
      ZNeu(4,2) = NMix_ZNeu(4,2)
      ZNeu(4,3) = NMix_ZNeu(4,3)
      ZNeu(4,4) = NMix_ZNeu(4,4)

* transforming the convention to complex mixing matrices
	
!	if (MNeu(1)*MNeu(2)*MNeu(3)*MNeu(4).lt.0d0) then
!	  do i=1,4
!		if (MNeu(i).lt.0d0) then
!		   MNeu(i) = -1d0*MNeu(i)
!		   do j=1,4
!			ZNeu(i,j) = cI*ZNeu(i,j)
!		   enddo	
!		endif 
!	  enddo
!	endif 

! ===== Chargino masses and mixing =====

      MCha(1) = Mass_MCha(1) 
      MCha(2) = Mass_MCha(2) 

      do i=1,2
            MCha2(i) = MCha(i)**2 
      enddo

      UCha(1,1) = UMix_UCha(1,1)
      UCha(1,2) = UMix_UCha(1,2)
      UCha(2,1) = UMix_UCha(2,1)
      UCha(2,2) = UMix_UCha(2,2)

      VCha(1,1) = VMix_VCha(1,1)
      VCha(1,2) = VMix_VCha(1,2)
      VCha(2,1) = VMix_VCha(2,1)
      VCha(2,2) = VMix_VCha(2,2)

* transforming the convention to complex mixing matrices

      if (MCha(1)*MCha(2).lt.0d0) then

      do i=1,2
      if (MCha(i).lt.0d0) then
      MCha(i) = -1d0*MCha(i)
      do j=1,2
      UCha(i,j) = cI*UCha(i,j)
      VCha(i,j) = cI*VCha(i,j)
      enddo	
      endif 
      enddo
      endif

      endif

! ===== Neutralino and Chargino mass together =====

      MNeuCha(1) = Mass_MNeu(1) 
      MNeuCha(2) = Mass_MNeu(2) 
      MNeuCha(3) = Mass_MNeu(3) 
      MNeuCha(4) = Mass_MNeu(4) 
      MNeuCha(5) = Mass_MCha(1)
      MNeuCha(6) = Mass_MCha(2)
	
	
! ===== Higgs masses and widths ====
      if(rerror.eq.1) call DMNLO_HiggsMasses

      MHk0(1) = Mh0
      MHk0(2) = MHH
      MHk0(3) = MA0
      MHk0(4) = MG0

      wHk0(1) = wh0
      wHk0(2) = wHH
      wHk0(3) = wA0
      wHk0(4) = wG0

      MHk(1) = Mh0
      MHk(2) = MHH
      MHk(3) = MA0
      MHk(4) = MG0
      MHk(5) = MHp
      MHk(6) = MGp
            
      wHk(1) = wh0
      wHk(2) = wHH
      wHk(3) = wA0
      wHk(4) = wG0
      wHk(5) = wHp
      wHk(6) = wGp

	end



! =================================================================
      subroutine DMNLO_SfermionMasses
      implicit none

#include "DMNLO_Model.h"


      integer t, g
      double precision i3part, qpart

      external NoMixing, Eigen2x2

      i3part = .5D0*MZ2*C2B
      qpart = SW2*MZ2*C2B
      DSf(1,1) = +i3part
      DSf(1,2) = 0
      DSf(2,2) = Qe*qpart
      DSf(1,2) = -i3part - DSf(2,2)
      DSf(2,3) = Qu*qpart
      DSf(1,3) = +i3part - DSf(2,3)
      DSf(2,4) = Qd*qpart
      DSf(1,4) = -i3part - DSf(2,4)

	do g = 1, 3
      MSS(1,2,g) = MSL(g)
      MSS(2,2,g) = MSE(g)
      MSS(1,3,g) = MSQ(g)
      MSS(2,3,g) = MSU(g)
      MSS(1,4,g) = MSQ(g)
      MSS(2,4,g) = MSD(g)

      MSS2(1,2,g) = MSS(1,2,g)**2
      MSS2(2,2,g) = MSS(2,2,g)**2
      MSS2(1,3,g) = MSS(1,3,g)**2
      MSS2(2,3,g) = MSS(2,3,g)**2
      MSS2(1,4,g) = MSS(1,4,g)**2
      MSS2(2,4,g) = MSS(2,4,g)**2

      Xf(2,g) = Af(2,g) - dconjg(MUE)*TB
      Xf(3,g) = Af(3,g) - dconjg(MUE)/TB
      Xf(4,g) = Af(4,g) - dconjg(MUE)*TB

      call NoMixing(MSf(1,1,g), MSf2(1,1,g), USf(1,1,1,g),
     &      MSS2(1,2,g) + DSf(1,1), 0D0, dcmplx(0D0), *999)

      do t = 2, 4
      call Eigen2x2(MSf(1,t,g), MSf2(1,t,g), USf(1,1,t,g),
     &        MSS2(1,t,g) + DSf(1,t) + Mf2(t,g),
     &        MSS2(2,t,g) + DSf(2,t) + Mf2(t,g),
     &        Mf(t,g)*dconjg(Xf(t,g)), *999)
      enddo
      enddo


999	return 
	end


! =================================================================
	subroutine DMNLO_CharginoMasses
	implicit none

#include "DMNLO_Model.h"

	double complex X(2,2)
	external SingularValues

	X(1,1) = M_2
	X(2,1) = sqrt2*MW*CB
	X(1,2) = sqrt2*MW*SB
	X(2,2) = MUE

	call SingularValues(2, 2, X, 2, MCha, UCha, 2, VCha, 2)

	MCha2(1) = MCha(1)**2
	MCha2(2) = MCha(2)**2
	
	end


! =================================================================
	subroutine DMNLO_NeutralinoMasses
	implicit none

#include "DMNLO_Model.h"

	double complex Y(4,4)
	external PseudoEigensystem

	Y(1,1) = M_1
	Y(2,2) = M_2
	Y(3,1) = -MZ*SW*CB
	Y(1,3) = Y(3,1)
	Y(4,1) = MZ*SW*SB
	Y(1,4) = Y(4,1)
	Y(3,2) = MZ*CW*CB
	Y(2,3) = Y(3,2)
	Y(4,2) = -MZ*CW*SB
	Y(2,4) = Y(4,2)
	Y(4,3) = -MUE
	Y(3,4) = Y(4,3)
	Y(2,1) = 0
	Y(1,2) = 0
	Y(3,3) = 0
	Y(4,4) = 0

	call PseudoEigensystem(4, Y, 4, MNeu, ZNeu, 4)

	MNeu2(1) = MNeu(1)**2
	MNeu2(2) = MNeu(2)**2
	MNeu2(3) = MNeu(3)**2
	MNeu2(4) = MNeu(4)**2
	end

! =================================================================
	subroutine DMNLO_HiggsMasses
	implicit none

#include "DMNLO_Model.h"

	double precision avgMH2, deltaMH2, TA

	MA02 = MA0**2

	avgMH2 = .5D0*(MA02 + MZ2)
	deltaMH2 = sqrt(avgMH2**2 - MA02*MZ2*C2B**2)
	Mh02 = avgMH2 - deltaMH2
	MHH2 = avgMH2 + deltaMH2

	MHp2 = MA02 + MW2

	Mh0 = sqrt(Mh02)

	MHH = sqrt(MHH2)

	MHp = sqrt(MHp2)

	TA = (-(MA02 + MZ2)*SB*CB)/
     &    (MZ2*CB2 + MA02*SB2 - Mh02)
	CA2 = 1/(1 + TA**2)
	SA2 = 1 - CA2
	SA = -sqrt(SA2)
	CA = SA/TA

	S2A = 2*SA*CA
	C2A = CA2 - SA2
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA


	end

! =================================================================
	subroutine DMNLO_ModelDigest

	implicit none

#include "DMNLO_Model.h"

	print *, "============================"
	print *, "SM parameters"
	print *, "----------------------------"
	print *, "1/Alfa =", 1/Alfa
	print *, "Alfas =", Alfas
	print *, "SW2 =", SW2
	print *, "----------------------------"
	print *, "Gauge-boson masses:"
	print *, "MW =", MW,"  wW =",wW
	print *, "MZ =", MZ,"  wZ =",wZ
	print *, "----------------------------"
	print *, "Fermion masses:"
	print *, "ME =", ME
	print *, "MM =", MM
	print *, "ML =", ML
	print *, "MU =", MU
	print *, "MC =", MC
	print *, "MT =", MT
	print *, "MD =", MD
	print *, "MS =", MS
	print *, "MB =", MB
	print *, "============================"
	print *, "MSSM parameters"
	print *, "----------------------------"
	print *, "TB =", TB
	print *, "M_1 =", M_1
	print *, "M_2 =", M_2
	print *, "MUE =", MUE
	print *, "MSusy =", MSusy
	print *, "Atau =", Af(2,3)
	print *, "At =", Af(3,3)
	print *, "Ab =", Af(4,3)
	print *, "----------------------------"
	print *, "Higgs masses:"
	print *, "(using tree-level Higgs masses)"
	print *, "Mh0 =", Mh0,"  wh0 =", wh0
	print *, "MHH =", MHH,"  wHH =", wHH
	print *, "MA0 =", MA0,"  wA0 =", wA0
	print *, "MHp =", MHp
	print *, "SA =", SA
	print *, "CA =", CA
	print *, "----------------------------"
	print *, "Sfermion masses:"
	print *, "Sneutrino  e:  ", MSf(1,1,1)
	print *, "           mu: ", MSf(1,1,2)
	print *, "           tau:", MSf(1,1,3)
	print *, "Slepton    e:  ", MSf(1,2,1), MSf(2,2,1)
	print *, "           mu: ", MSf(1,2,2), MSf(2,2,2)
	print *, "           tau:", MSf(1,2,3), MSf(2,2,3)
	print *, "u-Squark   u:  ", MSf(1,3,1), MSf(2,3,1)
	print *, "           c:  ", MSf(1,3,2), MSf(2,3,2)
	print *, "           t:  ", MSf(1,3,3), MSf(2,3,3)
	print *, "d-Squark   d:  ", MSf(1,4,1), MSf(2,4,1)
	print *, "           s:  ", MSf(1,4,2), MSf(2,4,2)
	print *, "           b:  ", MSf(1,4,3), MSf(2,4,3)
	print *, "----------------------------"
	print *, "Chargino masses:"
	print *, "MCha(1) =", MCha(1)
	print *, "MCha(2) =", MCha(2)
	print *, "----------------------------"
	print *, "Neutralino masses:"
	print *, "MNeu(1) =", MNeu(1)
	print *, "MNeu(2) =", MNeu(2)
	print *, "MNeu(3) =", MNeu(3)
	print *, "MNeu(4) =", MNeu(4)
	print *, "----------------------------"
	print *, "Gluino mass:"
	print *, "MGl =", MGl
	print *, "============================"

	end

