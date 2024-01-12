		double complex AL(6,2,4,3),AR(6,2,4,3),B(2,4),C(2,4)
		double complex FL(8,4,3),FR(8,4,3),KL(6,2,4,3),KR(6,2,4,3),
     -                 GLR(8,2,2,4,3), G(8,2,2,4,3), HL(6,6,8), HR(6,6,8) 
    	double complex BL(6,4),BR(6,4),CL(2,2,4,3,4),DL(6,6,4,4),DR(6,6,4,4)
		double complex S(2,2,2,2,4,4,3)
		double complex Aij(2,2,4,3)
		double complex CLL(2,2,4,3,4)
		double complex JL(2,4,3),JR(2,4,3)

	
		common /mycoupl/ AL,AR,B,C
		common /mycoupl/ FL,FR,GLR,G,HL,HR,KL,KR
		common /mycoupl/ BL,BR,CL,DL,DR,S,Aij,CLL
		common /mycoupl/ JL,JR
