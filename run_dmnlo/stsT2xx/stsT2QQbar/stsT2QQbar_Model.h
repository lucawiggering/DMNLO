* model.h

#include "../../util/DMNLO_Model.h"

* here one can include process-specific parameters

      double precision zerocheck
      double precision corrProp, corrVert, corrBox
      double complex stst2qqbar_M2Treesingleterms(5,5)
      double complex stst2qqbar_M2Treesinglet,stst2qqbar_M2Treeoctet

      parameter (zerocheck = 1.D-30)

      double precision wNeu(4), wCha(2)

      common /add_para/ wNeu, wCha
      common /Coulomb/ corrProp, corrVert, corrBox
      common /Tree/ stst2qqbar_M2Treesingleterms
      common /colorsplit/ stst2qqbar_M2Treesinglet,stst2qqbar_M2Treeoctet

