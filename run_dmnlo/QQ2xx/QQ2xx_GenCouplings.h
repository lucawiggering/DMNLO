	double complex ApL,BpL,AL(6,2,4,4,3,3),AR(6,2,4,4,3,3),CL(6),CR(6),DL(6,4,4,3,3),DR(6,4,4,3,3),DLt(6,4,4,3,3),DRt(6,4,4,3,3),DLu(6,4,4,3,3),DRu(6,4,4,3,3),ELc(4,2,4,3),ER(4,2,4,3),FL(4,2,4,3),FR(4,2,4,3),
     -                 GL(4,2,4,3), GR(4,2,4,3), HL(4,2,4,3), HR(4,2,4,3),ILR(6,6,6),JLR(4,6,6),LLR(4,4,6),VVV(4,4,4),GGH(4,4,6),GGV(4,4,4),sF4a(2,2,2,2,4,4,3),
     -                 sF4b(2,2,2,2,4,4,3),sF4(2,2,2,2,4,4,3),GLR(6,2,2,4,3),KLR(6,6,2,2,4,3)
	double complex DLM(4), DRM(4),ELcM(2,4,3),ERM(2,4,3),FLM(2,4,3),FRM(2,4,3),
     -                 GLM(2,4,3), GRM(2,4,3), HLM(2,4,3), HRM(2,4,3),BL(4,4,4,3,3),BR(4,4,4,3,3),BLt(4,4,4,3,3),BRt(4,4,4,3,3),BLu(4,4,4,3,3),BRu(4,4,4,3,3),BV(4),BA(4),GluL(2,4,4,3,3),
     -                 GluR(2,4,4,3,3),KL(6,2,4,3),KR(6,2,4,3) 
	double complex DLho(4), DRho(4),ELcho(2,4,3),ERho(2,4,3),FLho(2,4,3),FRho(2,4,3),
     -                 GLho(2,4,3), GRho(2,4,3), HLho(2,4,3), HRho(2,4,3) 
	double complex HSfSf(6,2,2,4,4,3,3), VSfSf(4,2,2,4,4,3,3),HSfSft(6,2,2,4,4,3,3),HSfSfu(6,2,2,4,4,3,3),VSfSft(4,2,2,4,4,3,3),VSfSfu(4,2,2,4,4,3,3),HHSfSf(6,6,2,2,4,4,3),VVSfSf(4,4,2,2,4,4,3),
     -                 VVSfSft(4,4,2,2,4,4,3),VVSfSfu(4,4,2,2,4,4,3)
	double complex ELbox(4,2,4,3), ERbox(4,2,4,3), FLbox(4,2,4,3), FRbox(4,2,4,3),
     -                 GLbox(4,2,4,3), GRbox(4,2,4,3), HLbox(4,2,4,3), HRbox(4,2,4,3) 
	double precision MHv2(4)


	common /mycoupl/ ApL,BpL,AL,AR,CL,CR,DL,DR,DLt,DRt,DLu,DRu,ELc,ER,FL,FR,GL,GR,HL,HR,HSfSf,VSfSf,ILR,JLR,HSfSft,HSfSfu,VSfSft,VSfSfu,HHSfSf,VVSfSf,VVSfSft,VVSfSfu,LLR,VVV,GGH,GGV,sF4a,sF4b,sF4,GLR,KLR
	common /mycoupl/ DLM, DRM, ELcM, ERM, FLM, FRM, GLM, GRM, HLM, HRM, BL, BR, BLt, BRt, BLu, BRu, BA, BV, GluL, GluR, KL, KR
	common /mycoupl/ DLho, DRho, ELcho, ERho, FLho, FRho, GLho, GRho, HLho, HRho
	common /mymass/ MHv2
