* model.h

#include "../../util/DMNLO_Model.h"

* here one can include process-specific parameters

! The parameter FPSki turns on (FPSki=1) and off (FPSki=0d0) the full polarization sum for the gluon with momentum ki.
! In the case of FPSk1=FPSk2=0d0, the longitudinal polarizations have to be subtracted through ghosts.

      

        double precision zerocheck

        parameter (zerocheck = 1.D-30)

	double precision wNeu(4), wCha(2), FPSk1,FPSk2
	common /add_para/ wNeu, wCha,FPSk1,FPSk2
