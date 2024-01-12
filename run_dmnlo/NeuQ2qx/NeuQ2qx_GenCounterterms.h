#include "../util/DMNLO_RenScheme.h"

	double complex dAL(6,2,4,3),dAR(6,2,4,3),dFL(8,4,3),dFR(8,4,3),dKL(6,2,4,3),dKR(6,2,4,3)
	double complex dGLR(8,2,2,4,3), dG(8,2,2,4,3)
	double complex dBL(4,4,3), dBR(4,4,3), dCL(2,2,4,3,4)
        double complex Prop(2,2,4,3), PropVL(4,3), PropVR(4,3), PropSL(4,3), PropSR(4,3)

	common /countercoupl/ dAL, dAR, dFL, dFR, dKL, dKR, dGLR, dG
	common /countercoupl/ dBL, dBR, dCL
	common /countercoupl/ Prop,PropVL,PropVR,PropSL,PropSR
