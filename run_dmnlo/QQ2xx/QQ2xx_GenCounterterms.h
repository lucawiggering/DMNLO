#include "../util/DMNLO_RenScheme.h"

	double complex dAL(6,2,4,3),dAR(6,2,4,3),dFL(8,4,3),dFR(8,4,3),dKL(6,2,4,3),dKR(6,2,4,3),dKLR(6,6,2,2,4,3),dHHSfSf(6,6,2,2,4,4,3),dVVSfSf(4,4,2,2,4,4,3)
	double complex dGLR(8,2,2,4,3), dHSfSf(6,2,2,4,4,3,3),dHSfSft(6,2,2,4,4,3,3),dHSfSfu(6,2,2,4,4,3,3)
	double complex dBL(4,4,3), dBR(4,4,3), dVSfSf(4,2,2,4,4,3,3),dVSfSft(4,2,2,4,4,3,3),dVSfSfu(4,2,2,4,4,3,3)
        double complex Prop(2,2,4,3), PropVL(4,3), PropVR(4,3), PropSL(4,3), PropSR(4,3)

	common /countercoupl/ dAL, dAR, dFL, dFR, dKL, dKR, dGLR, dHSfSf, dHSfSft, dHSfSfu, dKLR, dHHSfSf, dVVSfSf
	common /countercoupl/ dBL,dBR,dVSfSf,dVSfSft,dVSfSfu
	common /countercoupl/ Prop,PropVL,PropVR,PropSL,PropSR
