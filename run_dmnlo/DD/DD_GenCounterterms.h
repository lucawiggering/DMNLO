#include "../util/DMNLO_RenScheme.h"

	double complex dBL(4,3),dBR(4,3),dBpL(4,3),dBpR(4,3),dDL(6,4,3),dDR(6,4,3),dELc(4,2,4,3),dER(4,2,4,3),dFL(4,2,4,3),dFR(4,2,4,3),
     -                 dGL(4,2,4,3), dGR(4,2,4,3), dHL(4,2,4,3), dHR(4,2,4,3),dEFT_SI_L(4,2,4,3),dEFT_SI_R(4,2,4,3),dEFT_SD_L(4,2,4,3),dEFT_SD_R(4,2,4,3)
	double complex DDPropCounter(2,2,4,3)
	double complex dPifVR_SUSY(4,3), dPifVL_SUSY(4,3), dPifSR_SUSY(4,3), dPifSL_SUSY(4,3)
	double complex PifVR_SUSY(4,3), PifVL_SUSY(4,3), PifSR_SUSY(4,3), PifSL_SUSY(4,3)
	double complex dZfL_SUSY(4,3), dZfR_SUSY(4,3)

	common /countercoupl/ dBL,dBR,dDL,dDR,dELc,dER,dFL,dFR,dGL,dGR,dHL,dHR,dEFT_SI_L,dEFT_SI_R,dEFT_SD_L,dEFT_SD_R
	common /counterprop/ DDPropCounter
	common /wavefunction/ PifVR_SUSY, PifVL_SUSY, PifSR_SUSY, PifSL_SUSY, dPifVR_SUSY, dPifVL_SUSY, dPifSR_SUSY, dPifSL_SUSY, dZfL_SUSY, dZfR_SUSY
