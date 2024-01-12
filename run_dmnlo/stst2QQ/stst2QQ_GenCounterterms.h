#include "../util/DMNLO_RenScheme.h"

	double complex dVSfSf(4,2,2,4,4,3,3), Prop(2,2,4,3), dHSfSf(6,2,2,4,4,3,3), dGLR(6,2,2,4,3)
	double complex dBL(4,4,4,3,3), dBR(4,4,4,3,3)
	double complex dGluL(2,4,4,3,3),dGluR(2,4,4,3,3),dGluLt1(2,4,4,3,3),dGluRt1(2,4,4,3,3),dGluLt2(2,4,4,3,3),dGluRt2(2,4,4,3,3)
	double complex dAL(6,2,4,4,3,3),dAR(6,2,4,4,3,3),dER(4,2,4,3),dELc(4,2,4,3)
	double complex dALt1(6,2,4,4,3,3),dARt1(6,2,4,4,3,3),dALt2(6,2,4,4,3,3),dARt2(6,2,4,4,3,3)
	double complex dDL(6,4,4,3,3), dDR(6,4,4,3,3)

	common /countercoupl/ dVSfSf,Prop,dHSfSf,dGLR
	common /countercoupl/ dBL,dBR
	common /countercoupl/ dGluL,dGluR,dGluLt1,dGluRt1,dGluLt2,dGluRt2
	common /countercoupl/ dAL,dAR,dER,dELc
	common /countercoupl/ dALt1,dARt1,dALt2,dARt2
	common /countercoupl/ dDL, dDR
