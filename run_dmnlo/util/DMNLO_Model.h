* Model file for all processes 
* (DO NOT write here process specific stuff used only in one code!)

* model.h

	double precision pi, sqrt2, zeta3, hbar_c2

	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (sqrt2 = 1.41421356237309504880168872421D0)
	parameter (zeta3 = 1.2020569031595942D0)
	
* hbar c^2 in picobarn
	 parameter (hbar_c2 = 3.8937966D8)


	double complex cI
	parameter (cI = (0D0, 1D0))

	double precision Qe, Qu, Qd
	parameter (Qe = -1, Qu = 2/3D0, Qd = -1/3D0)

	double precision NC, CF, CV, TF
	parameter (NC = 3)
	parameter (CF = 4/3D0)
	parameter (CV = 3D0)
	parameter (TF = 1/2D0)

* DRbar Scale
	double precision QScale, Qscalealphas
	double precision Divergence
	double precision Esoftmax

	common /renorm/ QScale, Qscalealphas, Divergence, Esoftmax

* other switches
	
	integer rerror
      integer imass, ischeme
	
	common /iswitch/ rerror, imass, ischeme
	
* SM parameters

	double precision MZ, MZ2, MW, MW2, SW, S2W, C2W, TW, TW2, CW, CW2, SW2, wZ, wW
	double precision Alfa, Alfa2

	double precision ME, ME2, MM, MM2, ML, ML2
	double precision MU, MU2, MC, MC2, MT, MT2
	double precision MD, MD2, MS, MS2, MB, MB2
	double precision MBSM, MBrun, MBos, mbQdr
	double precision MTrun, MTos
	double precision MTscale,MBscale
	
	double complex CKM(3,3)
	double precision Mf(4,3), Mf2(4,3), Mfos(4,3), MfH(4,3,7), Mfos2(4,3), Mfm(5), Mfho(5), dMb
	double precision MG0, MG02, MGp, MGp2
	double precision wG0, wGp, width(4)
	double precision EL, GSMZ, AlfasMZ, Alfas2MZ, GS, Alfas, Alfas2
	double precision g1h1, g2h1, g1h2, g2h2

	common /sm_para/ CKM
	common /sm_para/ MTscale,MBscale
	common /sm_para/ Mf, Mf2, Mfos, MfH, Mfos2, Mfm, Mfho , dMb
	common /sm_para/ ME, ME2, MM, MM2, ML, ML2
	common /sm_para/ MU, MU2, MC, MC2, MT, MT2
	common /sm_para/ MD, MD2, MS, MS2, MB, MB2
	common /sm_para/ MBSM, MBrun, MBos, MTrun, MTos, mbQdr
	common /sm_para/ MG0, MG02, MGp, MGp2, wG0, wGp, width
	common /sm_para/ Alfa, Alfa2, EL, GS, Alfas, Alfas2, GSMZ, AlfasMZ, Alfas2MZ
	common /sm_para/ MZ, MZ2, MW, MW2, CW, CW2, SW, SW2, S2W, C2W, TW, TW2, wZ, wW
	common /sm_para/ g1h1, g2h1, g1h2, g2h2


* MSSM parameters

	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex USf(2,2,4,3),USfmix(2,2,4,3), Af(2:4,3), Xf(2:4,3)
	double complex Atau, At, Ab, MUE
	double precision MCha(2), MCha2(2), MNeu(4), MNeu2(4), MNeuCha(6)
	double precision MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	double precision MSE(3), MSL(3), MSU(3), MSD(3), MSQ(3)
	double precision MSE2(3), MSL2(3), MSU2(3), MSD2(3), MSQ2(3)
	double precision MSf(2,4,3), WSf(2,4,3), MSf2(2,4,3), MSusy
	double precision Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	double precision wh0, wHH, wA0, wHp
	double precision wt, wst(2)
	double precision M_1, M_2, MGl, MGl2
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA

	common /mssm_para/ UCha, VCha, ZNeu
	common /mssm_para/ USf, USfmix, Af, Xf
	common /mssm_para/ Atau, At, Ab, MUE
	common /mssm_para/ MCha, MCha2, MNeu, MNeu2, MNeuCha
	common /mssm_para/ MSS, MSS2, DSf
	common /mssm_para/ MSE, MSL, MSU, MSD, MSQ
	common /mssm_para/ MSE2, MSL2, MSU2, MSD2, MSQ2
	common /mssm_para/ MSf, WSf, MSf2, MSusy
	common /mssm_para/ Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2, wh0, wHH, wA0, wHp
	common /mssm_para/ wt, wst
	common /mssm_para/ M_1, M_2, MGl, MGl2
	common /mssm_para/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssm_para/ CA, SA, CA2, SA2, C2A, S2A
	common /mssm_para/ CAB, SAB, CBA, SBA

* shortcuts for 2->3 process

	double precision MHk0(4), wHk0(4), MHk(6), wHk(6)
	
	common /para_short/ MHk0, wHk0, MHk, wHk
	
#define CKMC(a,b) dconjg(CKM(a,b))
#define USfC(a,b,t,g) dconjg(USf(a,b,t,g))
#define VChaC(a,b) dconjg(VCha(a,b))
#define UChaC(a,b) dconjg(UCha(a,b))
#define ZNeuC(a,b) dconjg(ZNeu(a,b))



