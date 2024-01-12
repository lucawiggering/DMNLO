	double complex BL(4,4,4,3,3), BR(4,4,4,3,3), BLt(4,4,4,3,3), BRt(4,4,4,3,3), BLu(4,4,4,3,3), BRu(4,4,4,3,3),BV(4),BA(4)
	double complex DL(6,4,4,3,3), DR(6,4,4,3,3), DLt(6,4,4,3,3), DRt(6,4,4,3,3),DLu(6,4,4,3,3), DRu(6,4,4,3,3)


	double complex GLR(6,2,2,4,3), HSfSf(6,2,2,4,4,3,3),HSf(2)
	double complex VSfSf(4,2,2,4,4,3,3), VSfSft(4,2,2,4,4,3,3), VSfSfu(4,2,2,4,4,3,3),C(2,4)
	double complex GluL(2,4,4,3,3), GluR(2,4,4,3,3),hff(2)
	double complex GluLt1(2,4,4,3,3), GluRt1(2,4,4,3,3)
	double complex GluLt2(2,4,4,3,3), GluRt2(2,4,4,3,3)

	common /mycoupl/ BL, BR, BRu, BLu, BRt, BLt
	common /mycoupl/ DL, DR, DLt, DRt, DLu, DRu


	common /mycoupl/ GLR,C,HSf,hff
	common /mycoupl/ VSfSf, VSfSft, VSfSfu
	common /mycoupl/ GluL, GluR
	common /mycoupl/ GluLt1, GluRt1
	common /mycoupl/ GluLt2, GluRt2
