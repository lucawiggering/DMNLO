	double complex BL(4,4,4,3,3), BR(4,4,4,3,3), BLt(4,4,4,3,3), BRt(4,4,4,3,3), BLu(4,4,4,3,3), BRu(4,4,4,3,3),BV(4),BA(4)
	double complex DL(6,4,4,3,3), DR(6,4,4,3,3), DLt(6,4,4,3,3), DRt(6,4,4,3,3),DLu(6,4,4,3,3), DRu(6,4,4,3,3)
	double complex AL(6,2,4,4,3,3), AR(6,2,4,4,3,3)
	double complex ALt1(6,2,4,4,3,3), ARt1(6,2,4,4,3,3)
	double complex ALt2(6,2,4,4,3,3), ARt2(6,2,4,4,3,3)
	double complex ALt1ho(6,2,4,4,3,3), ARt1ho(6,2,4,4,3,3)
	double complex GluL(2,4,4,3,3), GluR(2,4,4,3,3)
	double complex GluLt1(2,4,4,3,3), GluRt1(2,4,4,3,3)
	double complex GluLt2(2,4,4,3,3), GluRt2(2,4,4,3,3)
	double complex C(2,4), VSfSf(4,2,2,4,4,3,3), VSfSft(4,2,2,4,4,3,3), VSfSfu(4,2,2,4,4,3,3),SfSfgg(2,4,3,2,4,3),SfSfVV(4,2,4,3,2,4,3)
	double complex GLR(6,2,2,4,3), HSfSf(6,2,2,4,4,3,3), HSfSft(4,2,2,4,4,3,3), HSfSfu(4,2,2,4,4,3,3)
	double complex Aij(2,2,4,4,3,3)
	



	common /mycoupl/ BL, BR, BRu, BLu, BRt, BLt, BV, BA
	common /mycoupl/ DL, DR, DLt, DRt, DLu, DRu
	common /mycoupl/ AL, AR
	common /mycoupl/ ALt1, ARt1
	common /mycoupl/ ALt2, ARt2
	common /mycoupl/ ALt1ho, ARt1ho
	common /mycoupl/ GluL, GluR
	common /mycoupl/ GluLt1, GluRt1
	common /mycoupl/ GluLt2, GluRt2
	common /mycoupl/ C, VSfSf, VSfSft, VSfSfu,SfSfgg,SfSfVV
	common /mycoupl/ GLR, HSfSf, HSfSft, HSfSfu
	common /mycoupl/ Aij
	
