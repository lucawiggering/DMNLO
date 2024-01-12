	double complex ApL(4,3),BpL(4,3),AL(4,3),AR(4,3),BL(4,3),BR(4,3),CL(6,4,3),CR(6,4,3),DL(6,4,3),DR(6,4,3),
     -				ELc(4,2,4,3),ER(4,2,4,3),FL(4,2,4,3),FR(4,2,4,3),GL(4,2,4,3), GR(4,2,4,3), HL(4,2,4,3), HR(4,2,4,3) 
	double complex DLho(6,4,3), DRho(6,4,3), ELcho(4,2,4,3), ERho(4,2,4,3), FLho(4,2,4,3), FRho(4,2,4,3), GLho(4,2,4,3),
     -			    GRho(4,2,4,3), HLho(4,2,4,3), HRho(4,2,4,3) 
	double complex HSfSf(6,2,2,4,3), GSfSf(2,2,4,3), ZSfSf(2,2,4,3)
	double complex IL(2,4,3), IR(2,4,3), JL(2,4,3), JR(2,4,3)

	common /mycoupl/ ApL,BpL,AL,AR,BL,BR,CL,CR,DL,DR,ELc,ER,FL,FR,GL,GR,HL,HR,IL,IR,JL,JR,HSfSf,GSfSf,ZSfSf
	common /mycoupl/ DLho,DRho,ELcho,ERho,FLho,FRho,GLho,GRho,HLho,HRho
