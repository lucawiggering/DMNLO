	double complex ApL,BpL,AL,AR,BL,BR,CL(6),CR(6),DL(6),DR(6),ELc(4,2,4,3),ER(4,2,4,3),FL(4,2,4,3),FR(4,2,4,3),
     -                 GL(4,2,4,3), GR(4,2,4,3), HL(4,2,4,3), HR(4,2,4,3) 
	double complex DLM(4), DRM(4),ELcM(2,4,3),ERM(2,4,3),FLM(2,4,3),FRM(2,4,3),
     -                 GLM(2,4,3), GRM(2,4,3), HLM(2,4,3), HRM(2,4,3) 
	double complex DLho(6), DRho(6),ELcho(4,2,4,3),ERho(4,2,4,3),FLho(4,2,4,3),FRho(4,2,4,3),
     -                 GLho(4,2,4,3), GRho(4,2,4,3), HLho(4,2,4,3), HRho(4,2,4,3) 
	double complex HSfSf(6,2,2,4,3), GSfSf(2,2,4,3), ZSfSf(2,2,4,3), WSfSf(2,2,4,3)
	double precision MHv2(4)


	common /mycoupl/ ApL,BpL,AL,AR,BL,BR,CL,CR,DL,DR,ELc,ER,FL,FR,GL,GR,HL,HR,HSfSf,GSfSf,ZSfSf,WSfSf
	common /mycoupl/ DLM, DRM, ELcM, ERM, FLM, FRM, GLM, GRM, HLM, HRM
	common /mycoupl/ DLho, DRho, ELcho, ERho, FLho, FRho, GLho, GRho, HLho, HRho
	common /mymass/ MHv2
