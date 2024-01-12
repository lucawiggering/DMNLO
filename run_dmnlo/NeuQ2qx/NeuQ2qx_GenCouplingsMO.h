		double complex AL_MO(6,2,4,3),AR_MO(6,2,4,3),B_MO(2,4),C_MO(2,4)
		double complex FL_MO(8,4,3),FR_MO(8,4,3),KL_MO(6,2,4,3),KR_MO(6,2,4,3)
                double complex GLR_MO(8,2,2,4,3), G_MO(8,2,2,4,3), HL_MO(6,6,8), HR_MO(6,6,8) 
    	        double complex BL_MO(6,4),BR_MO(6,4),CL_MO(2,2,4,3,4),DL_MO(6,6,4,4),DR_MO(6,6,4,4)
		double complex S_MO(2,2,2,2,4,4,3)
		double precision CLL_MO(2,2,4,3,4)
		double precision JL_MO(2,4,3),JR_MO(2,4,3)

	
		common /mycouplMO/ AL_MO,AR_MO,B_MO,C_MO
		common /mycouplMO/ FL_MO,FR_MO,GLR_MO,G_MO,HL_MO,HR_MO,KL_MO,KR_MO
		common /mycouplMO/ BL_MO,BR_MO,CL_MO,DL_MO,DR_MO,S_MO,CLL_MO
		common /mycouplMO/ JL_MO,JR_MO
