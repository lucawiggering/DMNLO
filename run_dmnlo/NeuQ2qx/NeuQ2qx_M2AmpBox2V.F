! --- MsBoxV  ---	

	double precision function MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)        
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc
        double precision mq1c,xv,wquark
              
        MsBoxV = 0d0
        MsBoxV = dreal((0.0015831434944115277*(A1Lc*BLc*mq1c*(-8*Bs0L*m1**2 + 4*Bs1R*m1**3 + 4*Bs2R*m1*m2**2 - 8*Bv0R*m1*m2**2 + 4*Bv1L*m1**2*m2**2 - 8*Bv2L*m1**2*m2**2 - 4*Bv2L*m2**4 + 
     &     32*Bs0R*m1*m3 + 8*Bs1L*m1**2*m3 - 4*Bs2L*m1**2*m3 + 4*Bs3L*m1**2*m3 - 16*Bv0L*m1**2*m3 - 4*Bv1R*m1**3*m3 - 4*Bs2L*m2**2*m3 - 16*Bv0L*m2**2*m3 - 
     &     4*Bv1R*m1*m2**2*m3 + 8*Bv2R*m1*m2**2*m3 + 4*Bv3R*m1*m2**2*m3 - 8*Bs0L*m3**2 + 4*Bs1R*m1*m3**2 + 4*Bs2R*m1*m3**2 + 8*Bs3R*m1*m3**2 - 
     &     8*Bv0R*m1*m3**2 + 4*Bv1L*m1**2*m3**2 - 4*Bv2L*m1**2*m3**2 - 4*Bv3L*m1**2*m3**2 - 8*Bv2L*m2**2*m3**2 - 4*Bv3L*m2**2*m3**2 + 4*Bs3L*m3**3 + 
     &     4*Bv3R*m1*m3**3 + 4*Bv2L*m2**2*sman + 4*Bs2L*m3*sman + 16*Bv0L*m3*sman + 4*Bv1R*m1*m3*sman + 4*Bv2L*m3**2*sman + 4*Bv3L*m3**2*sman + 8*Bs0L*tman - 
     &     4*Bs1R*m1*tman + 4*Bv2L*m2**2*tman - 4*Bs3L*m3*tman - 4*Bs2R*m1*uman + 8*Bv0R*m1*uman - 4*Bv1L*m1**2*uman + 4*Bv2L*m1**2*uman + 4*Bv2L*m2**2*uman - 
     &     4*Bv3R*m1*m3*uman - 4*Bv2L*sman*uman - (2*Bv1R*m1**3*m2**2*m3 + 2*Bv2R*m1*m2**4*m3 - 4*Bs0L*m1**2*m3**2 - 2*Bs1R*m1**3*m3**2 + 2*Bv1L*m1**4*m3**2 - 
     &        2*Bs2R*m1*m2**2*m3**2 - 4*Bv0R*m1*m2**2*m3**2 + 2*Bv2L*m1**2*m2**2*m3**2 - 2*Bs3L*m1**2*m3**3 - 2*Bv3R*m1*m2**2*m3**3 + 2*Bs3R*m1*m3**4 - 
     &        2*Bv3L*m1**2*m3**4 - 8*Bs0L*m1**2*m4**2 - 2*Bs1R*m1**3*m4**2 - 2*Bs2R*m1*m2**2*m4**2 - 8*Bv0R*m1*m2**2*m4**2 + 2*Bv1L*m1**2*m2**2*m4**2 + 
     &        2*Bv2L*m2**4*m4**2 + 8*Bs0R*m1*m3*m4**2 - 2*Bs3L*m1**2*m3*m4**2 - 4*Bv0L*m1**2*m3*m4**2 + 2*Bv1R*m1**3*m3*m4**2 - 4*Bv0L*m2**2*m3*m4**2 + 
     &        2*Bv1R*m1*m2**2*m3*m4**2 + 4*Bv2R*m1*m2**2*m3*m4**2 - 2*Bv3R*m1*m2**2*m3*m4**2 - 8*Bs0L*m3**2*m4**2 - 2*Bs1R*m1*m3**2*m4**2 - 
     &        2*Bs2R*m1*m3**2*m4**2 + 4*Bs3R*m1*m3**2*m4**2 - 8*Bv0R*m1*m3**2*m4**2 + 2*Bv1L*m1**2*m3**2*m4**2 + 2*Bv2L*m1**2*m3**2*m4**2 - 
     &        2*Bv3L*m1**2*m3**2*m4**2 - 2*Bv3L*m2**2*m3**2*m4**2 - 2*Bs3L*m3**3*m4**2 - 2*Bv3R*m1*m3**3*m4**2 - 4*Bs0L*m4**4 - 2*Bs1R*m1*m4**4 - 
     &        2*Bs2R*m1*m4**4 + 2*Bs3R*m1*m4**4 - 4*Bv0R*m1*m4**4 + 2*Bv1L*m2**2*m4**4 + 2*Bv2L*m2**2*m4**4 - 2*Bv3L*m2**2*m4**4 - 2*Bs3L*m3*m4**4 + 
     &        2*Bv1R*m1*m3*m4**4 + 2*Bv2R*m1*m3*m4**4 - 2*Bv3R*m1*m3*m4**4 + 4*Bs0L*m1**2*sman + 2*Bs1R*m1**3*sman - Bv1L*m1**4*sman + 2*Bs2R*m1*m2**2*sman + 
     &        4*Bv0R*m1*m2**2*sman - Bv1L*m1**2*m2**2*sman - Bv2L*m1**2*m2**2*sman - Bv2L*m2**4*sman + 2*Bs3L*m1**2*m3*sman + 2*Bv3R*m1*m2**2*m3*sman - 
     &        4*Bs3R*m1*m3**2*sman - Bv1L*m1**2*m3**2*sman + 3*Bv3L*m1**2*m3**2*sman - Bv2L*m2**2*m3**2*sman + Bv3L*m2**2*m3**2*sman + Bv3L*m3**4*sman + 
     &        4*Bs0L*m4**2*sman + 2*Bs1R*m1*m4**2*sman + 2*Bs2R*m1*m4**2*sman - 4*Bs3R*m1*m4**2*sman + 4*Bv0R*m1*m4**2*sman - 2*Bv1L*m1**2*m4**2*sman - 
     &        Bv2L*m1**2*m4**2*sman + Bv3L*m1**2*m4**2*sman - Bv1L*m2**2*m4**2*sman - 2*Bv2L*m2**2*m4**2*sman + 3*Bv3L*m2**2*m4**2*sman + 
     &        2*Bs3L*m3*m4**2*sman + 4*Bv0L*m3*m4**2*sman + 2*Bv3R*m1*m3*m4**2*sman - Bv1L*m3**2*m4**2*sman - Bv2L*m3**2*m4**2*sman + 
     &        2*Bv3L*m3**2*m4**2*sman - Bv1L*m4**4*sman - Bv2L*m4**4*sman + Bv3L*m4**4*sman + 2*Bs3R*m1*sman**2 + Bv1L*m1**2*sman**2 - Bv3L*m1**2*sman**2 + 
     &        Bv2L*m2**2*sman**2 - Bv3L*m2**2*sman**2 - 2*Bv3L*m3**2*sman**2 + Bv1L*m4**2*sman**2 + Bv2L*m4**2*sman**2 - 2*Bv3L*m4**2*sman**2 + Bv3L*sman**3 + 
     &        Bv1L*m1**4*tman + Bv1L*m1**2*m2**2*tman + Bv2L*m1**2*m2**2*tman + Bv2L*m2**4*tman - 2*Bv1R*m1**3*m3*tman - 4*Bv2R*m1*m2**2*m3*tman + 
     &        2*Bs2R*m1*m3**2*tman + 4*Bv0R*m1*m3**2*tman + Bv1L*m1**2*m3**2*tman - 2*Bv2L*m1**2*m3**2*tman - Bv3L*m1**2*m3**2*tman + Bv2L*m2**2*m3**2*tman - 
     &        Bv3L*m2**2*m3**2*tman + 2*Bv3R*m1*m3**3*tman - Bv3L*m3**4*tman + 4*Bs0L*m4**2*tman + 2*Bs2R*m1*m4**2*tman + 4*Bv0R*m1*m4**2*tman + 
     &        2*Bv1L*m1**2*m4**2*tman + Bv2L*m1**2*m4**2*tman - Bv3L*m1**2*m4**2*tman + Bv1L*m2**2*m4**2*tman - Bv3L*m2**2*m4**2*tman - 
     &        2*Bv1R*m1*m3*m4**2*tman - 4*Bv2R*m1*m3*m4**2*tman + 2*Bv3R*m1*m3*m4**2*tman + Bv1L*m3**2*m4**2*tman + Bv2L*m3**2*m4**2*tman - 
     &        2*Bv3L*m3**2*m4**2*tman + Bv1L*m4**4*tman + Bv2L*m4**4*tman - Bv3L*m4**4*tman - 2*Bs2R*m1*sman*tman - 4*Bv0R*m1*sman*tman + 
     &        Bv2L*m1**2*sman*tman + Bv3L*m1**2*sman*tman + Bv2L*m2**2*sman*tman + Bv3L*m2**2*sman*tman - 2*Bv3R*m1*m3*sman*tman + Bv2L*m3**2*sman*tman + 
     &        Bv3L*m3**2*sman*tman + Bv2L*m4**2*sman*tman + Bv3L*m4**2*sman*tman - Bv2L*sman**2*tman - Bv1L*m1**2*tman**2 - Bv2L*m1**2*tman**2 - 
     &        2*Bv2L*m2**2*tman**2 + 2*Bv2R*m1*m3*tman**2 - Bv2L*m3**2*tman**2 + Bv3L*m3**2*tman**2 - Bv1L*m4**2*tman**2 - 2*Bv2L*m4**2*tman**2 + 
     &        Bv3L*m4**2*tman**2 - Bv3L*sman*tman**2 + Bv2L*tman**3 + 2*Bs2L*m3*(m2**2 + m4**2 - tman)*(m1**2 + m4**2 - uman) + 
     &        2*Bs1L*m3*(m1**2 + m4**2 - uman)**2 - (2*Bv1R*m1*m2**2*m3 - 4*Bs0L*m3**2 - 2*Bs1R*m1*m3**2 - Bv3L*m1**2*m3**2 - Bv3L*m2**2*m3**2 - 
     &           2*Bs3L*m3**3 - Bv3L*m3**4 - 4*Bs0L*m4**2 - 2*Bs1R*m1*m4**2 - 4*Bv0R*m1*m4**2 - Bv3L*m1**2*m4**2 - Bv3L*m2**2*m4**2 - 2*Bs3L*m3*m4**2 + 
     &           2*Bv1R*m1*m3*m4**2 - 2*Bv3L*m3**2*m4**2 - Bv3L*m4**4 + 4*Bs0L*sman + 2*Bs1R*m1*sman + Bv3L*m1**2*sman + Bv3L*m2**2*sman + 2*Bs3L*m3*sman + 
     &           Bv3L*m3**2*sman + Bv3L*m4**2*sman + Bv2L*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1R*m1*m3*tman + 
     &           Bv1L*(m1**4 + m3**2*m4**2 + m4**4 - m3**2*sman - m4**2*sman + sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + 3*m3**2 + 2*m4**2 - sman + tman) + m2**2*(3*m4**2 - sman + tman)))*uman + 
     &        (Bv1L*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - Bv3L*(m3**2 + m4**2 - sman) + Bv2L*(m2**2 + m4**2 - tman))*uman**2 - Bv1L*uman**3)*xv) + 
     &  A1Rc*BRc*mq1c*(-8*Bs0R*m1**2 + 4*Bs1L*m1**3 + 4*Bs2L*m1*m2**2 - 8*Bv0L*m1*m2**2 + 4*Bv1R*m1**2*m2**2 - 8*Bv2R*m1**2*m2**2 - 4*Bv2R*m2**4 + 
     &     32*Bs0L*m1*m3 + 8*Bs1R*m1**2*m3 - 4*Bs2R*m1**2*m3 + 4*Bs3R*m1**2*m3 - 16*Bv0R*m1**2*m3 - 4*Bv1L*m1**3*m3 - 4*Bs2R*m2**2*m3 - 16*Bv0R*m2**2*m3 - 
     &     4*Bv1L*m1*m2**2*m3 + 8*Bv2L*m1*m2**2*m3 + 4*Bv3L*m1*m2**2*m3 - 8*Bs0R*m3**2 + 4*Bs1L*m1*m3**2 + 4*Bs2L*m1*m3**2 + 8*Bs3L*m1*m3**2 - 
     &     8*Bv0L*m1*m3**2 + 4*Bv1R*m1**2*m3**2 - 4*Bv2R*m1**2*m3**2 - 4*Bv3R*m1**2*m3**2 - 8*Bv2R*m2**2*m3**2 - 4*Bv3R*m2**2*m3**2 + 4*Bs3R*m3**3 + 
     &     4*Bv3L*m1*m3**3 + 4*Bv2R*m2**2*sman + 4*Bs2R*m3*sman + 16*Bv0R*m3*sman + 4*Bv1L*m1*m3*sman + 4*Bv2R*m3**2*sman + 4*Bv3R*m3**2*sman + 8*Bs0R*tman - 
     &     4*Bs1L*m1*tman + 4*Bv2R*m2**2*tman - 4*Bs3R*m3*tman - 4*Bs2L*m1*uman + 8*Bv0L*m1*uman - 4*Bv1R*m1**2*uman + 4*Bv2R*m1**2*uman + 4*Bv2R*m2**2*uman - 
     &     4*Bv3L*m1*m3*uman - 4*Bv2R*sman*uman - (2*Bv1L*m1**3*m2**2*m3 + 2*Bv2L*m1*m2**4*m3 - 4*Bs0R*m1**2*m3**2 - 2*Bs1L*m1**3*m3**2 + 2*Bv1R*m1**4*m3**2 - 
     &        2*Bs2L*m1*m2**2*m3**2 - 4*Bv0L*m1*m2**2*m3**2 + 2*Bv2R*m1**2*m2**2*m3**2 - 2*Bs3R*m1**2*m3**3 - 2*Bv3L*m1*m2**2*m3**3 + 2*Bs3L*m1*m3**4 - 
     &        2*Bv3R*m1**2*m3**4 - 8*Bs0R*m1**2*m4**2 - 2*Bs1L*m1**3*m4**2 - 2*Bs2L*m1*m2**2*m4**2 - 8*Bv0L*m1*m2**2*m4**2 + 2*Bv1R*m1**2*m2**2*m4**2 + 
     &        2*Bv2R*m2**4*m4**2 + 8*Bs0L*m1*m3*m4**2 - 2*Bs3R*m1**2*m3*m4**2 - 4*Bv0R*m1**2*m3*m4**2 + 2*Bv1L*m1**3*m3*m4**2 - 4*Bv0R*m2**2*m3*m4**2 + 
     &        2*Bv1L*m1*m2**2*m3*m4**2 + 4*Bv2L*m1*m2**2*m3*m4**2 - 2*Bv3L*m1*m2**2*m3*m4**2 - 8*Bs0R*m3**2*m4**2 - 2*Bs1L*m1*m3**2*m4**2 - 
     &        2*Bs2L*m1*m3**2*m4**2 + 4*Bs3L*m1*m3**2*m4**2 - 8*Bv0L*m1*m3**2*m4**2 + 2*Bv1R*m1**2*m3**2*m4**2 + 2*Bv2R*m1**2*m3**2*m4**2 - 
     &        2*Bv3R*m1**2*m3**2*m4**2 - 2*Bv3R*m2**2*m3**2*m4**2 - 2*Bs3R*m3**3*m4**2 - 2*Bv3L*m1*m3**3*m4**2 - 4*Bs0R*m4**4 - 2*Bs1L*m1*m4**4 - 
     &        2*Bs2L*m1*m4**4 + 2*Bs3L*m1*m4**4 - 4*Bv0L*m1*m4**4 + 2*Bv1R*m2**2*m4**4 + 2*Bv2R*m2**2*m4**4 - 2*Bv3R*m2**2*m4**4 - 2*Bs3R*m3*m4**4 + 
     &        2*Bv1L*m1*m3*m4**4 + 2*Bv2L*m1*m3*m4**4 - 2*Bv3L*m1*m3*m4**4 + 4*Bs0R*m1**2*sman + 2*Bs1L*m1**3*sman - Bv1R*m1**4*sman + 2*Bs2L*m1*m2**2*sman + 
     &        4*Bv0L*m1*m2**2*sman - Bv1R*m1**2*m2**2*sman - Bv2R*m1**2*m2**2*sman - Bv2R*m2**4*sman + 2*Bs3R*m1**2*m3*sman + 2*Bv3L*m1*m2**2*m3*sman - 
     &        4*Bs3L*m1*m3**2*sman - Bv1R*m1**2*m3**2*sman + 3*Bv3R*m1**2*m3**2*sman - Bv2R*m2**2*m3**2*sman + Bv3R*m2**2*m3**2*sman + Bv3R*m3**4*sman + 
     &        4*Bs0R*m4**2*sman + 2*Bs1L*m1*m4**2*sman + 2*Bs2L*m1*m4**2*sman - 4*Bs3L*m1*m4**2*sman + 4*Bv0L*m1*m4**2*sman - 2*Bv1R*m1**2*m4**2*sman - 
     &        Bv2R*m1**2*m4**2*sman + Bv3R*m1**2*m4**2*sman - Bv1R*m2**2*m4**2*sman - 2*Bv2R*m2**2*m4**2*sman + 3*Bv3R*m2**2*m4**2*sman + 
     &        2*Bs3R*m3*m4**2*sman + 4*Bv0R*m3*m4**2*sman + 2*Bv3L*m1*m3*m4**2*sman - Bv1R*m3**2*m4**2*sman - Bv2R*m3**2*m4**2*sman + 
     &        2*Bv3R*m3**2*m4**2*sman - Bv1R*m4**4*sman - Bv2R*m4**4*sman + Bv3R*m4**4*sman + 2*Bs3L*m1*sman**2 + Bv1R*m1**2*sman**2 - Bv3R*m1**2*sman**2 + 
     &        Bv2R*m2**2*sman**2 - Bv3R*m2**2*sman**2 - 2*Bv3R*m3**2*sman**2 + Bv1R*m4**2*sman**2 + Bv2R*m4**2*sman**2 - 2*Bv3R*m4**2*sman**2 + Bv3R*sman**3 + 
     &        Bv1R*m1**4*tman + Bv1R*m1**2*m2**2*tman + Bv2R*m1**2*m2**2*tman + Bv2R*m2**4*tman - 2*Bv1L*m1**3*m3*tman - 4*Bv2L*m1*m2**2*m3*tman + 
     &        2*Bs2L*m1*m3**2*tman + 4*Bv0L*m1*m3**2*tman + Bv1R*m1**2*m3**2*tman - 2*Bv2R*m1**2*m3**2*tman - Bv3R*m1**2*m3**2*tman + Bv2R*m2**2*m3**2*tman - 
     &        Bv3R*m2**2*m3**2*tman + 2*Bv3L*m1*m3**3*tman - Bv3R*m3**4*tman + 4*Bs0R*m4**2*tman + 2*Bs2L*m1*m4**2*tman + 4*Bv0L*m1*m4**2*tman + 
     &        2*Bv1R*m1**2*m4**2*tman + Bv2R*m1**2*m4**2*tman - Bv3R*m1**2*m4**2*tman + Bv1R*m2**2*m4**2*tman - Bv3R*m2**2*m4**2*tman - 
     &        2*Bv1L*m1*m3*m4**2*tman - 4*Bv2L*m1*m3*m4**2*tman + 2*Bv3L*m1*m3*m4**2*tman + Bv1R*m3**2*m4**2*tman + Bv2R*m3**2*m4**2*tman - 
     &        2*Bv3R*m3**2*m4**2*tman + Bv1R*m4**4*tman + Bv2R*m4**4*tman - Bv3R*m4**4*tman - 2*Bs2L*m1*sman*tman - 4*Bv0L*m1*sman*tman + 
     &        Bv2R*m1**2*sman*tman + Bv3R*m1**2*sman*tman + Bv2R*m2**2*sman*tman + Bv3R*m2**2*sman*tman - 2*Bv3L*m1*m3*sman*tman + Bv2R*m3**2*sman*tman + 
     &        Bv3R*m3**2*sman*tman + Bv2R*m4**2*sman*tman + Bv3R*m4**2*sman*tman - Bv2R*sman**2*tman - Bv1R*m1**2*tman**2 - Bv2R*m1**2*tman**2 - 
     &        2*Bv2R*m2**2*tman**2 + 2*Bv2L*m1*m3*tman**2 - Bv2R*m3**2*tman**2 + Bv3R*m3**2*tman**2 - Bv1R*m4**2*tman**2 - 2*Bv2R*m4**2*tman**2 + 
     &        Bv3R*m4**2*tman**2 - Bv3R*sman*tman**2 + Bv2R*tman**3 + 2*Bs2R*m3*(m2**2 + m4**2 - tman)*(m1**2 + m4**2 - uman) + 
     &        2*Bs1R*m3*(m1**2 + m4**2 - uman)**2 - (2*Bv1L*m1*m2**2*m3 - 4*Bs0R*m3**2 - 2*Bs1L*m1*m3**2 - Bv3R*m1**2*m3**2 - Bv3R*m2**2*m3**2 - 
     &           2*Bs3R*m3**3 - Bv3R*m3**4 - 4*Bs0R*m4**2 - 2*Bs1L*m1*m4**2 - 4*Bv0L*m1*m4**2 - Bv3R*m1**2*m4**2 - Bv3R*m2**2*m4**2 - 2*Bs3R*m3*m4**2 + 
     &           2*Bv1L*m1*m3*m4**2 - 2*Bv3R*m3**2*m4**2 - Bv3R*m4**4 + 4*Bs0R*sman + 2*Bs1L*m1*sman + Bv3R*m1**2*sman + Bv3R*m2**2*sman + 2*Bs3R*m3*sman + 
     &           Bv3R*m3**2*sman + Bv3R*m4**2*sman + Bv2R*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1L*m1*m3*tman + 
     &           Bv1R*(m1**4 + m3**2*m4**2 + m4**4 - m3**2*sman - m4**2*sman + sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + 3*m3**2 + 2*m4**2 - sman + tman) + m2**2*(3*m4**2 - sman + tman)))*uman + 
     &        (Bv1R*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - Bv3R*(m3**2 + m4**2 - sman) + Bv2R*(m2**2 + m4**2 - tman))*uman**2 - Bv1R*uman**3)*xv) + 
     &  A1Lc*BRc*(8*Bs1L*m1**4 + 8*Bs2L*m1**2*m2**2 + 8*Bv1R*m1**3*m2**2 + 8*Bv2R*m1*m2**4 + 32*Bs0L*m1**2*m3 + 8*Bs1R*m1**3*m3 - 8*Bv1L*m1**4*m3 + 
     &     8*Bs2R*m1*m2**2*m3 + 32*Bv0R*m1*m2**2*m3 - 16*Bv1L*m1**2*m2**2*m3 + 8*Bv2L*m1**2*m2**2*m3 - 8*Bs0R*m1*m3**2 + 8*Bs1L*m1**2*m3**2 + 
     &     8*Bs2L*m1**2*m3**2 + 8*Bs3L*m1**2*m3**2 + 8*Bv1R*m1**3*m3**2 + 8*Bs2L*m2**2*m3**2 + 8*Bv0L*m2**2*m3**2 + 8*Bv1R*m1*m2**2*m3**2 + 
     &     8*Bv2R*m1*m2**2*m3**2 + 8*Bv3R*m1*m2**2*m3**2 + 16*Bs0L*m3**3 + 4*Bs1R*m1*m3**3 + 4*Bs2R*m1*m3**3 + 4*Bs3R*m1*m3**3 + 16*Bv0R*m1*m3**3 - 
     &     8*Bv1L*m1**2*m3**3 - 4*Bv1L*m2**2*m3**3 + 4*Bv2L*m2**2*m3**3 - 4*Bv3L*m2**2*m3**3 + 4*Bs1L*m3**4 + 4*Bs2L*m3**4 + 4*Bs3L*m3**4 + 4*Bv1R*m1*m3**4 + 
     &     4*Bv2R*m1*m3**4 + 4*Bv3R*m1*m3**4 + 8*Bs0R*m1*m4**2 + 8*Bs1L*m1**2*m4**2 - 8*Bv0L*m2**2*m4**2 + 8*Bv2R*m1*m2**2*m4**2 + 16*Bs0L*m3*m4**2 + 
     &     4*Bs1R*m1*m3*m4**2 + 4*Bs2R*m1*m3*m4**2 - 4*Bs3R*m1*m3*m4**2 + 16*Bv0R*m1*m3*m4**2 - 8*Bv1L*m1**2*m3*m4**2 - 4*Bv1L*m2**2*m3*m4**2 + 
     &     4*Bv2L*m2**2*m3*m4**2 + 4*Bv3L*m2**2*m3*m4**2 + 4*Bs1L*m3**2*m4**2 + 4*Bs2L*m3**2*m4**2 + 4*Bs3L*m3**2*m4**2 + 4*Bv1R*m1*m3**2*m4**2 + 
     &     4*Bv2R*m1*m3**2*m4**2 + 4*Bv3R*m1*m3**2*m4**2 - 8*Bs0R*m1*sman - 4*Bs1L*m1**2*sman + 2*Bs2L*m1**2*sman + 4*Bv0L*m1**2*sman + 2*Bv1R*m1**3*sman + 
     &     2*Bs2L*m2**2*sman + 4*Bv0L*m2**2*sman + 2*Bv1R*m1*m2**2*sman - 4*Bv2R*m1*m2**2*sman + 4*Bs3R*m1*m3*sman + 8*Bv1L*m1**2*m3*sman - 
     &     2*Bv3L*m1**2*m3*sman - 2*Bv3L*m2**2*m3*sman - 2*Bs2L*m3**2*sman - 4*Bv0L*m3**2*sman - 2*Bv1R*m1*m3**2*sman + 4*Bv1L*m3**3*sman + 
     &     2*Bv3L*m3**3*sman + 2*Bs2L*m4**2*sman + 4*Bv0L*m4**2*sman + 2*Bv1R*m1*m4**2*sman + 4*Bv1L*m3*m4**2*sman - 2*Bv3L*m3*m4**2*sman - 2*Bs2L*sman**2 - 
     &     4*Bv0L*sman**2 - 2*Bv1R*m1*sman**2 + 2*Bv3L*m3*sman**2 - 12*Bs1L*m1**2*tman - 2*Bs2L*m1**2*tman - 4*Bv0L*m1**2*tman - 2*Bv1R*m1**3*tman - 
     &     6*Bs2L*m2**2*tman - 4*Bv0L*m2**2*tman - 6*Bv1R*m1*m2**2*tman - 4*Bv2R*m1*m2**2*tman - 16*Bs0L*m3*tman - 4*Bs1R*m1*m3*tman - 4*Bs2R*m1*m3*tman - 
     &     16*Bv0R*m1*m3*tman + 8*Bv1L*m1**2*m3*tman + 2*Bv3L*m1**2*m3*tman + 4*Bv1L*m2**2*m3*tman - 4*Bv2L*m2**2*m3*tman + 2*Bv3L*m2**2*m3*tman - 
     &     8*Bs1L*m3**2*tman - 6*Bs2L*m3**2*tman - 4*Bs3L*m3**2*tman - 4*Bv0L*m3**2*tman - 6*Bv1R*m1*m3**2*tman - 4*Bv2R*m1*m3**2*tman - 
     &     4*Bv3R*m1*m3**2*tman + 2*Bv3L*m3**3*tman - 4*Bs1L*m4**2*tman - 2*Bs2L*m4**2*tman - 4*Bv0L*m4**2*tman - 2*Bv1R*m1*m4**2*tman + 
     &     2*Bv3L*m3*m4**2*tman - 4*Bv1L*m3*sman*tman + 4*Bs1L*tman**2 + 2*Bs2L*tman**2 + 4*Bv0L*tman**2 + 2*Bv1R*m1*tman**2 - 2*Bv3L*m3*tman**2 - 
     &     4*Bs1L*m1**2*uman - 6*Bs2L*m1**2*uman + 4*Bv0L*m1**2*uman - 6*Bv1R*m1**3*uman - 2*Bs2L*m2**2*uman + 4*Bv0L*m2**2*uman - 2*Bv1R*m1*m2**2*uman - 
     &     12*Bv2R*m1*m2**2*uman - 16*Bs0L*m3*uman - 4*Bs1R*m1*m3*uman - 4*Bs2R*m1*m3*uman - 16*Bv0R*m1*m3*uman + 8*Bv1L*m1**2*m3*uman - 
     &     2*Bv3L*m1**2*m3*uman + 4*Bv1L*m2**2*m3*uman - 4*Bv2L*m2**2*m3*uman - 2*Bv3L*m2**2*m3*uman - 4*Bs1L*m3**2*uman - 6*Bs2L*m3**2*uman - 
     &     4*Bs3L*m3**2*uman + 4*Bv0L*m3**2*uman - 6*Bv1R*m1*m3**2*uman - 8*Bv2R*m1*m3**2*uman - 4*Bv3R*m1*m3**2*uman - 2*Bv3L*m3**3*uman - 
     &     2*Bs2L*m4**2*uman + 4*Bv0L*m4**2*uman - 2*Bv1R*m1*m4**2*uman - 4*Bv2R*m1*m4**2*uman - 2*Bv3L*m3*m4**2*uman - 4*Bv1L*m3*sman*uman + 
     &     4*Bs1L*tman*uman + 4*Bs2L*tman*uman + 4*Bv1R*m1*tman*uman + 4*Bv2R*m1*tman*uman + 2*Bs2L*uman**2 - 4*Bv0L*uman**2 + 2*Bv1R*m1*uman**2 + 
     &     4*Bv2R*m1*uman**2 + 2*Bv3L*m3*uman**2 + (4*Bv1R*m1**3*m2**2*m3**2 + 4*Bv2R*m1*m2**4*m3**2 + 2*Bs1R*m1**3*m3**3 + 2*Bs2R*m1*m2**2*m3**3 - 
     &        2*Bv1L*m1**2*m2**2*m3**3 - 2*Bv2L*m2**4*m3**3 - 4*Bs0R*m1*m3**4 - 4*Bs3L*m1**2*m3**4 + 2*Bv1R*m1**3*m3**4 + 4*Bv0L*m2**2*m3**4 + 
     &        2*Bv2R*m1*m2**2*m3**4 - 4*Bv3R*m1*m2**2*m3**4 - 2*Bs3R*m1*m3**5 + 2*Bv3L*m2**2*m3**5 - 2*Bs3L*m3**6 - 2*Bv3R*m1*m3**6 - 8*Bs0L*m1**2*m3*m4**2 - 
     &        2*Bs1R*m1**3*m3*m4**2 - 2*Bs2R*m1*m2**2*m3*m4**2 - 8*Bv0R*m1*m2**2*m3*m4**2 + 2*Bv1L*m1**2*m2**2*m3*m4**2 + 2*Bv2L*m2**4*m3*m4**2 + 
     &        4*Bs0R*m1*m3**2*m4**2 - 4*Bs3L*m1**2*m3**2*m4**2 + 2*Bv1R*m1**3*m3**2*m4**2 - 4*Bv0L*m2**2*m3**2*m4**2 + 4*Bv1R*m1*m2**2*m3**2*m4**2 + 
     &        6*Bv2R*m1*m2**2*m3**2*m4**2 - 4*Bv3R*m1*m2**2*m3**2*m4**2 - 4*Bs0L*m3**3*m4**2 + 2*Bs1R*m1*m3**3*m4**2 + 2*Bs2R*m1*m3**3*m4**2 - 
     &        4*Bv0R*m1*m3**3*m4**2 - 2*Bv1L*m2**2*m3**3*m4**2 - 2*Bv2L*m2**2*m3**3*m4**2 - 4*Bs3L*m3**4*m4**2 + 2*Bv1R*m1*m3**4*m4**2 + 
     &        2*Bv2R*m1*m3**4*m4**2 - 4*Bv3R*m1*m3**4*m4**2 - 4*Bs0L*m3*m4**4 - 2*Bs1R*m1*m3*m4**4 - 2*Bs2R*m1*m3*m4**4 + 2*Bs3R*m1*m3*m4**4 - 
     &        4*Bv0R*m1*m3*m4**4 + 2*Bv1L*m2**2*m3*m4**4 + 2*Bv2L*m2**2*m3*m4**4 - 2*Bv3L*m2**2*m3*m4**4 - 2*Bs3L*m3**2*m4**4 + 2*Bv1R*m1*m3**2*m4**4 + 
     &        2*Bv2R*m1*m3**2*m4**4 - 2*Bv3R*m1*m3**2*m4**4 - 2*Bv1R*m1**3*m2**2*sman - 2*Bv2R*m1*m2**4*sman - 2*Bs1R*m1**3*m3*sman + Bv1L*m1**4*m3*sman - 
     &        2*Bs2R*m1*m2**2*m3*sman + Bv1L*m1**2*m2**2*m3*sman + Bv2L*m1**2*m2**2*m3*sman + Bv2L*m2**4*m3*sman + 8*Bs0R*m1*m3**2*sman + 
     &        6*Bs3L*m1**2*m3**2*sman - 2*Bv0L*m1**2*m3**2*sman - 2*Bv1R*m1**3*m3**2*sman - 6*Bv0L*m2**2*m3**2*sman - 2*Bv2R*m1*m2**2*m3**2*sman + 
     &        6*Bv3R*m1*m2**2*m3**2*sman + 4*Bs3R*m1*m3**3*sman + Bv1L*m1**2*m3**3*sman - Bv3L*m1**2*m3**3*sman + Bv2L*m2**2*m3**3*sman - 
     &        3*Bv3L*m2**2*m3**3*sman + 4*Bs3L*m3**4*sman - 2*Bv0L*m3**4*sman + 4*Bv3R*m1*m3**4*sman - Bv3L*m3**5*sman + 4*Bs0R*m1*m4**2*sman + 
     &        2*Bs3L*m1**2*m4**2*sman - 4*Bv0L*m1**2*m4**2*sman - 2*Bv1R*m1*m2**2*m4**2*sman - 2*Bv2R*m1*m2**2*m4**2*sman + 2*Bv3R*m1*m2**2*m4**2*sman - 
     &        2*Bs1R*m1*m3*m4**2*sman - 2*Bs2R*m1*m3*m4**2*sman + Bv2L*m1**2*m3*m4**2*sman - Bv3L*m1**2*m3*m4**2*sman + Bv1L*m2**2*m3*m4**2*sman + 
     &        Bv3L*m2**2*m3*m4**2*sman + 4*Bs3L*m3**2*m4**2*sman + 2*Bv0L*m3**2*m4**2*sman - 2*Bv1R*m1*m3**2*m4**2*sman - 2*Bv2R*m1*m3**2*m4**2*sman + 
     &        4*Bv3R*m1*m3**2*m4**2*sman + Bv1L*m3**3*m4**2*sman + Bv2L*m3**3*m4**2*sman - Bv1L*m3*m4**4*sman - Bv2L*m3*m4**4*sman + Bv3L*m3*m4**4*sman - 
     &        4*Bs0R*m1*sman**2 - 2*Bs3L*m1**2*sman**2 + 2*Bv0L*m1**2*sman**2 + 2*Bv0L*m2**2*sman**2 - 2*Bv3R*m1*m2**2*sman**2 - 2*Bs3R*m1*m3*sman**2 - 
     &        Bv1L*m1**2*m3*sman**2 + Bv3L*m1**2*m3*sman**2 - Bv2L*m2**2*m3*sman**2 + Bv3L*m2**2*m3*sman**2 - 2*Bs3L*m3**2*sman**2 + 4*Bv0L*m3**2*sman**2 - 
     &        2*Bv3R*m1*m3**2*sman**2 + 2*Bv3L*m3**3*sman**2 + 2*Bv0L*m4**2*sman**2 - Bv1L*m3*m4**2*sman**2 - Bv2L*m3*m4**2*sman**2 - 2*Bv0L*sman**3 - 
     &        Bv3L*m3*sman**3 + Bv1L*m1**4*m3*tman + Bv1L*m1**2*m2**2*m3*tman + Bv2L*m1**2*m2**2*m3*tman + Bv2L*m2**4*m3*tman - 2*Bv0L*m1**2*m3**2*tman - 
     &        2*Bv1R*m1**3*m3**2*tman - 2*Bv0L*m2**2*m3**2*tman - 6*Bv2R*m1*m2**2*m3**2*tman - 2*Bs2R*m1*m3**3*tman + Bv1L*m1**2*m3**3*tman - 
     &        Bv3L*m1**2*m3**3*tman + 3*Bv2L*m2**2*m3**3*tman - Bv3L*m2**2*m3**3*tman + 2*Bs3L*m3**4*tman - 2*Bv0L*m3**4*tman - 2*Bv2R*m1*m3**4*tman + 
     &        2*Bv3R*m1*m3**4*tman - Bv3L*m3**5*tman + 4*Bs0L*m3*m4**2*tman + 2*Bs2R*m1*m3*m4**2*tman + 4*Bv0R*m1*m3*m4**2*tman + 2*Bv1L*m1**2*m3*m4**2*tman + 
     &        Bv2L*m1**2*m3*m4**2*tman - Bv3L*m1**2*m3*m4**2*tman + Bv1L*m2**2*m3*m4**2*tman - Bv3L*m2**2*m3*m4**2*tman + 2*Bs3L*m3**2*m4**2*tman - 
     &        2*Bv0L*m3**2*m4**2*tman - 2*Bv1R*m1*m3**2*m4**2*tman - 4*Bv2R*m1*m3**2*m4**2*tman + 2*Bv3R*m1*m3**2*m4**2*tman + Bv1L*m3**3*m4**2*tman + 
     &        Bv2L*m3**3*m4**2*tman - 2*Bv3L*m3**3*m4**2*tman + Bv1L*m3*m4**4*tman + Bv2L*m3*m4**4*tman - Bv3L*m3*m4**4*tman + 2*Bv0L*m1**2*sman*tman + 
     &        2*Bv0L*m2**2*sman*tman + 2*Bv2R*m1*m2**2*sman*tman + 2*Bs2R*m1*m3*sman*tman - Bv2L*m1**2*m3*sman*tman + Bv3L*m1**2*m3*sman*tman - 
     &        Bv2L*m2**2*m3*sman*tman + Bv3L*m2**2*m3*sman*tman - 4*Bs3L*m3**2*sman*tman + 2*Bv0L*m3**2*sman*tman + 2*Bv2R*m1*m3**2*sman*tman - 
     &        2*Bv3R*m1*m3**2*sman*tman - Bv2L*m3**3*sman*tman + Bv3L*m3**3*sman*tman - 2*Bs3L*m4**2*sman*tman + 2*Bv0L*m4**2*sman*tman + 
     &        Bv2L*m3*m4**2*sman*tman + Bv3L*m3*m4**2*sman*tman + 2*Bs3L*sman**2*tman + Bv2L*m3*sman**2*tman - Bv1L*m1**2*m3*tman**2 - Bv2L*m1**2*m3*tman**2 - 
     &        2*Bv2L*m2**2*m3*tman**2 + 2*Bv0L*m3**2*tman**2 + 2*Bv2R*m1*m3**2*tman**2 - Bv2L*m3**3*tman**2 + Bv3L*m3**3*tman**2 - Bv1L*m3*m4**2*tman**2 - 
     &        2*Bv2L*m3*m4**2*tman**2 + Bv3L*m3*m4**2*tman**2 - 2*Bv0L*sman*tman**2 - Bv3L*m3*sman*tman**2 + Bv2L*m3*tman**3 + 
     &        2*Bs2L*(m2**2 + m4**2 - tman)*(m3**4 + m1**2*(2*m3**2 - sman) + sman*tman + m3**2*(m4**2 - sman - tman - uman)) + 
     &        2*Bs1L*(m3**4 + m1**2*(2*m3**2 - sman) + sman*tman + m3**2*(m4**2 - sman - tman - uman))*(m1**2 + m4**2 - uman) - 
     &        (-2*Bv0L*m1**2*m3**2 + 2*Bv1R*m1**3*m3**2 - 2*Bv0L*m2**2*m3**2 + 4*Bv1R*m1*m2**2*m3**2 + 2*Bv2R*m1*m2**2*m3**2 + 2*Bs1R*m1*m3**3 - 
     &           Bv3L*m1**2*m3**3 - Bv3L*m2**2*m3**3 - 2*Bs3L*m3**4 - 2*Bv0L*m3**4 + 2*Bv1R*m1*m3**4 - 2*Bv3R*m1*m3**4 - Bv3L*m3**5 - 4*Bs0L*m3*m4**2 - 
     &           2*Bs1R*m1*m3*m4**2 - 4*Bv0R*m1*m3*m4**2 - Bv3L*m1**2*m3*m4**2 - Bv3L*m2**2*m3*m4**2 - 2*Bs3L*m3**2*m4**2 - 2*Bv0L*m3**2*m4**2 + 
     &           4*Bv1R*m1*m3**2*m4**2 + 2*Bv2R*m1*m3**2*m4**2 - 2*Bv3R*m1*m3**2*m4**2 - 2*Bv3L*m3**3*m4**2 - Bv3L*m3*m4**4 + 2*Bv0L*m1**2*sman - 
     &           2*Bv1R*m1**3*sman + 2*Bv0L*m2**2*sman - 2*Bv1R*m1*m2**2*sman - 2*Bv2R*m1*m2**2*sman - 2*Bs1R*m1*m3*sman + Bv3L*m1**2*m3*sman + 
     &           Bv3L*m2**2*m3*sman + 2*Bs3L*m3**2*sman + 2*Bv0L*m3**2*sman - 2*Bv1R*m1*m3**2*sman + 4*Bv3R*m1*m3**2*sman + Bv3L*m3**3*sman + 
     &           2*Bv0L*m4**2*sman - 2*Bv1R*m1*m4**2*sman - 2*Bv2R*m1*m4**2*sman + 2*Bv3R*m1*m4**2*sman + Bv3L*m3*m4**2*sman - 2*Bv3R*m1*sman**2 + 
     &           Bv2L*m3*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1R*m1*m3**2*tman - 2*Bv2R*m1*m3**2*tman + 2*Bv2R*m1*sman*tman + 
     &           Bv1L*m3*(m1**4 + m3**2*m4**2 + m4**4 + m3**2*sman - m4**2*sman - sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + m3**2 + 2*m4**2 + sman + tman) + m2**2*(-2*m3**2 + 3*m4**2 + sman + tman)))*uman + 
     &        (Bv1L*m3*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - m3*(m3*(2*Bv0L - 2*Bv1R*m1 + Bv3L*m3) + Bv3L*m4**2) + (2*Bv0L - 2*Bv1R*m1 + Bv3L*m3)*sman + 
     &           Bv2L*m3*(m2**2 + m4**2 - tman))*uman**2 - Bv1L*m3*uman**3)*xv) + 
     &  A1Rc*BLc*(8*Bs1R*m1**4 + 8*Bs2R*m1**2*m2**2 + 8*Bv1L*m1**3*m2**2 + 8*Bv2L*m1*m2**4 + 32*Bs0R*m1**2*m3 + 8*Bs1L*m1**3*m3 - 8*Bv1R*m1**4*m3 + 
     &     8*Bs2L*m1*m2**2*m3 + 32*Bv0L*m1*m2**2*m3 - 16*Bv1R*m1**2*m2**2*m3 + 8*Bv2R*m1**2*m2**2*m3 - 8*Bs0L*m1*m3**2 + 8*Bs1R*m1**2*m3**2 + 
     &     8*Bs2R*m1**2*m3**2 + 8*Bs3R*m1**2*m3**2 + 8*Bv1L*m1**3*m3**2 + 8*Bs2R*m2**2*m3**2 + 8*Bv0R*m2**2*m3**2 + 8*Bv1L*m1*m2**2*m3**2 + 
     &     8*Bv2L*m1*m2**2*m3**2 + 8*Bv3L*m1*m2**2*m3**2 + 16*Bs0R*m3**3 + 4*Bs1L*m1*m3**3 + 4*Bs2L*m1*m3**3 + 4*Bs3L*m1*m3**3 + 16*Bv0L*m1*m3**3 - 
     &     8*Bv1R*m1**2*m3**3 - 4*Bv1R*m2**2*m3**3 + 4*Bv2R*m2**2*m3**3 - 4*Bv3R*m2**2*m3**3 + 4*Bs1R*m3**4 + 4*Bs2R*m3**4 + 4*Bs3R*m3**4 + 4*Bv1L*m1*m3**4 + 
     &     4*Bv2L*m1*m3**4 + 4*Bv3L*m1*m3**4 + 8*Bs0L*m1*m4**2 + 8*Bs1R*m1**2*m4**2 - 8*Bv0R*m2**2*m4**2 + 8*Bv2L*m1*m2**2*m4**2 + 16*Bs0R*m3*m4**2 + 
     &     4*Bs1L*m1*m3*m4**2 + 4*Bs2L*m1*m3*m4**2 - 4*Bs3L*m1*m3*m4**2 + 16*Bv0L*m1*m3*m4**2 - 8*Bv1R*m1**2*m3*m4**2 - 4*Bv1R*m2**2*m3*m4**2 + 
     &     4*Bv2R*m2**2*m3*m4**2 + 4*Bv3R*m2**2*m3*m4**2 + 4*Bs1R*m3**2*m4**2 + 4*Bs2R*m3**2*m4**2 + 4*Bs3R*m3**2*m4**2 + 4*Bv1L*m1*m3**2*m4**2 + 
     &     4*Bv2L*m1*m3**2*m4**2 + 4*Bv3L*m1*m3**2*m4**2 - 8*Bs0L*m1*sman - 4*Bs1R*m1**2*sman + 2*Bs2R*m1**2*sman + 4*Bv0R*m1**2*sman + 2*Bv1L*m1**3*sman + 
     &     2*Bs2R*m2**2*sman + 4*Bv0R*m2**2*sman + 2*Bv1L*m1*m2**2*sman - 4*Bv2L*m1*m2**2*sman + 4*Bs3L*m1*m3*sman + 8*Bv1R*m1**2*m3*sman - 
     &     2*Bv3R*m1**2*m3*sman - 2*Bv3R*m2**2*m3*sman - 2*Bs2R*m3**2*sman - 4*Bv0R*m3**2*sman - 2*Bv1L*m1*m3**2*sman + 4*Bv1R*m3**3*sman + 
     &     2*Bv3R*m3**3*sman + 2*Bs2R*m4**2*sman + 4*Bv0R*m4**2*sman + 2*Bv1L*m1*m4**2*sman + 4*Bv1R*m3*m4**2*sman - 2*Bv3R*m3*m4**2*sman - 2*Bs2R*sman**2 - 
     &     4*Bv0R*sman**2 - 2*Bv1L*m1*sman**2 + 2*Bv3R*m3*sman**2 - 12*Bs1R*m1**2*tman - 2*Bs2R*m1**2*tman - 4*Bv0R*m1**2*tman - 2*Bv1L*m1**3*tman - 
     &     6*Bs2R*m2**2*tman - 4*Bv0R*m2**2*tman - 6*Bv1L*m1*m2**2*tman - 4*Bv2L*m1*m2**2*tman - 16*Bs0R*m3*tman - 4*Bs1L*m1*m3*tman - 4*Bs2L*m1*m3*tman - 
     &     16*Bv0L*m1*m3*tman + 8*Bv1R*m1**2*m3*tman + 2*Bv3R*m1**2*m3*tman + 4*Bv1R*m2**2*m3*tman - 4*Bv2R*m2**2*m3*tman + 2*Bv3R*m2**2*m3*tman - 
     &     8*Bs1R*m3**2*tman - 6*Bs2R*m3**2*tman - 4*Bs3R*m3**2*tman - 4*Bv0R*m3**2*tman - 6*Bv1L*m1*m3**2*tman - 4*Bv2L*m1*m3**2*tman - 
     &     4*Bv3L*m1*m3**2*tman + 2*Bv3R*m3**3*tman - 4*Bs1R*m4**2*tman - 2*Bs2R*m4**2*tman - 4*Bv0R*m4**2*tman - 2*Bv1L*m1*m4**2*tman + 
     &     2*Bv3R*m3*m4**2*tman - 4*Bv1R*m3*sman*tman + 4*Bs1R*tman**2 + 2*Bs2R*tman**2 + 4*Bv0R*tman**2 + 2*Bv1L*m1*tman**2 - 2*Bv3R*m3*tman**2 - 
     &     4*Bs1R*m1**2*uman - 6*Bs2R*m1**2*uman + 4*Bv0R*m1**2*uman - 6*Bv1L*m1**3*uman - 2*Bs2R*m2**2*uman + 4*Bv0R*m2**2*uman - 2*Bv1L*m1*m2**2*uman - 
     &     12*Bv2L*m1*m2**2*uman - 16*Bs0R*m3*uman - 4*Bs1L*m1*m3*uman - 4*Bs2L*m1*m3*uman - 16*Bv0L*m1*m3*uman + 8*Bv1R*m1**2*m3*uman - 
     &     2*Bv3R*m1**2*m3*uman + 4*Bv1R*m2**2*m3*uman - 4*Bv2R*m2**2*m3*uman - 2*Bv3R*m2**2*m3*uman - 4*Bs1R*m3**2*uman - 6*Bs2R*m3**2*uman - 
     &     4*Bs3R*m3**2*uman + 4*Bv0R*m3**2*uman - 6*Bv1L*m1*m3**2*uman - 8*Bv2L*m1*m3**2*uman - 4*Bv3L*m1*m3**2*uman - 2*Bv3R*m3**3*uman - 
     &     2*Bs2R*m4**2*uman + 4*Bv0R*m4**2*uman - 2*Bv1L*m1*m4**2*uman - 4*Bv2L*m1*m4**2*uman - 2*Bv3R*m3*m4**2*uman - 4*Bv1R*m3*sman*uman + 
     &     4*Bs1R*tman*uman + 4*Bs2R*tman*uman + 4*Bv1L*m1*tman*uman + 4*Bv2L*m1*tman*uman + 2*Bs2R*uman**2 - 4*Bv0R*uman**2 + 2*Bv1L*m1*uman**2 + 
     &     4*Bv2L*m1*uman**2 + 2*Bv3R*m3*uman**2 + (4*Bv1L*m1**3*m2**2*m3**2 + 4*Bv2L*m1*m2**4*m3**2 + 2*Bs1L*m1**3*m3**3 + 2*Bs2L*m1*m2**2*m3**3 - 
     &        2*Bv1R*m1**2*m2**2*m3**3 - 2*Bv2R*m2**4*m3**3 - 4*Bs0L*m1*m3**4 - 4*Bs3R*m1**2*m3**4 + 2*Bv1L*m1**3*m3**4 + 4*Bv0R*m2**2*m3**4 + 
     &        2*Bv2L*m1*m2**2*m3**4 - 4*Bv3L*m1*m2**2*m3**4 - 2*Bs3L*m1*m3**5 + 2*Bv3R*m2**2*m3**5 - 2*Bs3R*m3**6 - 2*Bv3L*m1*m3**6 - 8*Bs0R*m1**2*m3*m4**2 - 
     &        2*Bs1L*m1**3*m3*m4**2 - 2*Bs2L*m1*m2**2*m3*m4**2 - 8*Bv0L*m1*m2**2*m3*m4**2 + 2*Bv1R*m1**2*m2**2*m3*m4**2 + 2*Bv2R*m2**4*m3*m4**2 + 
     &        4*Bs0L*m1*m3**2*m4**2 - 4*Bs3R*m1**2*m3**2*m4**2 + 2*Bv1L*m1**3*m3**2*m4**2 - 4*Bv0R*m2**2*m3**2*m4**2 + 4*Bv1L*m1*m2**2*m3**2*m4**2 + 
     &        6*Bv2L*m1*m2**2*m3**2*m4**2 - 4*Bv3L*m1*m2**2*m3**2*m4**2 - 4*Bs0R*m3**3*m4**2 + 2*Bs1L*m1*m3**3*m4**2 + 2*Bs2L*m1*m3**3*m4**2 - 
     &        4*Bv0L*m1*m3**3*m4**2 - 2*Bv1R*m2**2*m3**3*m4**2 - 2*Bv2R*m2**2*m3**3*m4**2 - 4*Bs3R*m3**4*m4**2 + 2*Bv1L*m1*m3**4*m4**2 + 
     &        2*Bv2L*m1*m3**4*m4**2 - 4*Bv3L*m1*m3**4*m4**2 - 4*Bs0R*m3*m4**4 - 2*Bs1L*m1*m3*m4**4 - 2*Bs2L*m1*m3*m4**4 + 2*Bs3L*m1*m3*m4**4 - 
     &        4*Bv0L*m1*m3*m4**4 + 2*Bv1R*m2**2*m3*m4**4 + 2*Bv2R*m2**2*m3*m4**4 - 2*Bv3R*m2**2*m3*m4**4 - 2*Bs3R*m3**2*m4**4 + 2*Bv1L*m1*m3**2*m4**4 + 
     &        2*Bv2L*m1*m3**2*m4**4 - 2*Bv3L*m1*m3**2*m4**4 - 2*Bv1L*m1**3*m2**2*sman - 2*Bv2L*m1*m2**4*sman - 2*Bs1L*m1**3*m3*sman + Bv1R*m1**4*m3*sman - 
     &        2*Bs2L*m1*m2**2*m3*sman + Bv1R*m1**2*m2**2*m3*sman + Bv2R*m1**2*m2**2*m3*sman + Bv2R*m2**4*m3*sman + 8*Bs0L*m1*m3**2*sman + 
     &        6*Bs3R*m1**2*m3**2*sman - 2*Bv0R*m1**2*m3**2*sman - 2*Bv1L*m1**3*m3**2*sman - 6*Bv0R*m2**2*m3**2*sman - 2*Bv2L*m1*m2**2*m3**2*sman + 
     &        6*Bv3L*m1*m2**2*m3**2*sman + 4*Bs3L*m1*m3**3*sman + Bv1R*m1**2*m3**3*sman - Bv3R*m1**2*m3**3*sman + Bv2R*m2**2*m3**3*sman - 
     &        3*Bv3R*m2**2*m3**3*sman + 4*Bs3R*m3**4*sman - 2*Bv0R*m3**4*sman + 4*Bv3L*m1*m3**4*sman - Bv3R*m3**5*sman + 4*Bs0L*m1*m4**2*sman + 
     &        2*Bs3R*m1**2*m4**2*sman - 4*Bv0R*m1**2*m4**2*sman - 2*Bv1L*m1*m2**2*m4**2*sman - 2*Bv2L*m1*m2**2*m4**2*sman + 2*Bv3L*m1*m2**2*m4**2*sman - 
     &        2*Bs1L*m1*m3*m4**2*sman - 2*Bs2L*m1*m3*m4**2*sman + Bv2R*m1**2*m3*m4**2*sman - Bv3R*m1**2*m3*m4**2*sman + Bv1R*m2**2*m3*m4**2*sman + 
     &        Bv3R*m2**2*m3*m4**2*sman + 4*Bs3R*m3**2*m4**2*sman + 2*Bv0R*m3**2*m4**2*sman - 2*Bv1L*m1*m3**2*m4**2*sman - 2*Bv2L*m1*m3**2*m4**2*sman + 
     &        4*Bv3L*m1*m3**2*m4**2*sman + Bv1R*m3**3*m4**2*sman + Bv2R*m3**3*m4**2*sman - Bv1R*m3*m4**4*sman - Bv2R*m3*m4**4*sman + Bv3R*m3*m4**4*sman - 
     &        4*Bs0L*m1*sman**2 - 2*Bs3R*m1**2*sman**2 + 2*Bv0R*m1**2*sman**2 + 2*Bv0R*m2**2*sman**2 - 2*Bv3L*m1*m2**2*sman**2 - 2*Bs3L*m1*m3*sman**2 - 
     &        Bv1R*m1**2*m3*sman**2 + Bv3R*m1**2*m3*sman**2 - Bv2R*m2**2*m3*sman**2 + Bv3R*m2**2*m3*sman**2 - 2*Bs3R*m3**2*sman**2 + 4*Bv0R*m3**2*sman**2 - 
     &        2*Bv3L*m1*m3**2*sman**2 + 2*Bv3R*m3**3*sman**2 + 2*Bv0R*m4**2*sman**2 - Bv1R*m3*m4**2*sman**2 - Bv2R*m3*m4**2*sman**2 - 2*Bv0R*sman**3 - 
     &        Bv3R*m3*sman**3 + Bv1R*m1**4*m3*tman + Bv1R*m1**2*m2**2*m3*tman + Bv2R*m1**2*m2**2*m3*tman + Bv2R*m2**4*m3*tman - 2*Bv0R*m1**2*m3**2*tman - 
     &        2*Bv1L*m1**3*m3**2*tman - 2*Bv0R*m2**2*m3**2*tman - 6*Bv2L*m1*m2**2*m3**2*tman - 2*Bs2L*m1*m3**3*tman + Bv1R*m1**2*m3**3*tman - 
     &        Bv3R*m1**2*m3**3*tman + 3*Bv2R*m2**2*m3**3*tman - Bv3R*m2**2*m3**3*tman + 2*Bs3R*m3**4*tman - 2*Bv0R*m3**4*tman - 2*Bv2L*m1*m3**4*tman + 
     &        2*Bv3L*m1*m3**4*tman - Bv3R*m3**5*tman + 4*Bs0R*m3*m4**2*tman + 2*Bs2L*m1*m3*m4**2*tman + 4*Bv0L*m1*m3*m4**2*tman + 2*Bv1R*m1**2*m3*m4**2*tman + 
     &        Bv2R*m1**2*m3*m4**2*tman - Bv3R*m1**2*m3*m4**2*tman + Bv1R*m2**2*m3*m4**2*tman - Bv3R*m2**2*m3*m4**2*tman + 2*Bs3R*m3**2*m4**2*tman - 
     &        2*Bv0R*m3**2*m4**2*tman - 2*Bv1L*m1*m3**2*m4**2*tman - 4*Bv2L*m1*m3**2*m4**2*tman + 2*Bv3L*m1*m3**2*m4**2*tman + Bv1R*m3**3*m4**2*tman + 
     &        Bv2R*m3**3*m4**2*tman - 2*Bv3R*m3**3*m4**2*tman + Bv1R*m3*m4**4*tman + Bv2R*m3*m4**4*tman - Bv3R*m3*m4**4*tman + 2*Bv0R*m1**2*sman*tman + 
     &        2*Bv0R*m2**2*sman*tman + 2*Bv2L*m1*m2**2*sman*tman + 2*Bs2L*m1*m3*sman*tman - Bv2R*m1**2*m3*sman*tman + Bv3R*m1**2*m3*sman*tman - 
     &        Bv2R*m2**2*m3*sman*tman + Bv3R*m2**2*m3*sman*tman - 4*Bs3R*m3**2*sman*tman + 2*Bv0R*m3**2*sman*tman + 2*Bv2L*m1*m3**2*sman*tman - 
     &        2*Bv3L*m1*m3**2*sman*tman - Bv2R*m3**3*sman*tman + Bv3R*m3**3*sman*tman - 2*Bs3R*m4**2*sman*tman + 2*Bv0R*m4**2*sman*tman + 
     &        Bv2R*m3*m4**2*sman*tman + Bv3R*m3*m4**2*sman*tman + 2*Bs3R*sman**2*tman + Bv2R*m3*sman**2*tman - Bv1R*m1**2*m3*tman**2 - Bv2R*m1**2*m3*tman**2 - 
     &        2*Bv2R*m2**2*m3*tman**2 + 2*Bv0R*m3**2*tman**2 + 2*Bv2L*m1*m3**2*tman**2 - Bv2R*m3**3*tman**2 + Bv3R*m3**3*tman**2 - Bv1R*m3*m4**2*tman**2 - 
     &        2*Bv2R*m3*m4**2*tman**2 + Bv3R*m3*m4**2*tman**2 - 2*Bv0R*sman*tman**2 - Bv3R*m3*sman*tman**2 + Bv2R*m3*tman**3 + 
     &        2*Bs2R*(m2**2 + m4**2 - tman)*(m3**4 + m1**2*(2*m3**2 - sman) + sman*tman + m3**2*(m4**2 - sman - tman - uman)) + 
     &        2*Bs1R*(m3**4 + m1**2*(2*m3**2 - sman) + sman*tman + m3**2*(m4**2 - sman - tman - uman))*(m1**2 + m4**2 - uman) - 
     &        (-2*Bv0R*m1**2*m3**2 + 2*Bv1L*m1**3*m3**2 - 2*Bv0R*m2**2*m3**2 + 4*Bv1L*m1*m2**2*m3**2 + 2*Bv2L*m1*m2**2*m3**2 + 2*Bs1L*m1*m3**3 - 
     &           Bv3R*m1**2*m3**3 - Bv3R*m2**2*m3**3 - 2*Bs3R*m3**4 - 2*Bv0R*m3**4 + 2*Bv1L*m1*m3**4 - 2*Bv3L*m1*m3**4 - Bv3R*m3**5 - 4*Bs0R*m3*m4**2 - 
     &           2*Bs1L*m1*m3*m4**2 - 4*Bv0L*m1*m3*m4**2 - Bv3R*m1**2*m3*m4**2 - Bv3R*m2**2*m3*m4**2 - 2*Bs3R*m3**2*m4**2 - 2*Bv0R*m3**2*m4**2 + 
     &           4*Bv1L*m1*m3**2*m4**2 + 2*Bv2L*m1*m3**2*m4**2 - 2*Bv3L*m1*m3**2*m4**2 - 2*Bv3R*m3**3*m4**2 - Bv3R*m3*m4**4 + 2*Bv0R*m1**2*sman - 
     &           2*Bv1L*m1**3*sman + 2*Bv0R*m2**2*sman - 2*Bv1L*m1*m2**2*sman - 2*Bv2L*m1*m2**2*sman - 2*Bs1L*m1*m3*sman + Bv3R*m1**2*m3*sman + 
     &           Bv3R*m2**2*m3*sman + 2*Bs3R*m3**2*sman + 2*Bv0R*m3**2*sman - 2*Bv1L*m1*m3**2*sman + 4*Bv3L*m1*m3**2*sman + Bv3R*m3**3*sman + 
     &           2*Bv0R*m4**2*sman - 2*Bv1L*m1*m4**2*sman - 2*Bv2L*m1*m4**2*sman + 2*Bv3L*m1*m4**2*sman + Bv3R*m3*m4**2*sman - 2*Bv3L*m1*sman**2 + 
     &           Bv2R*m3*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1L*m1*m3**2*tman - 2*Bv2L*m1*m3**2*tman + 2*Bv2L*m1*sman*tman + 
     &           Bv1R*m3*(m1**4 + m3**2*m4**2 + m4**4 + m3**2*sman - m4**2*sman - sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + m3**2 + 2*m4**2 + sman + tman) + m2**2*(-2*m3**2 + 3*m4**2 + sman + tman)))*uman + 
     &        (Bv1R*m3*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - m3*(m3*(2*Bv0R - 2*Bv1L*m1 + Bv3R*m3) + Bv3R*m4**2) + (2*Bv0R - 2*Bv1L*m1 + Bv3R*m3)*sman + 
     &           Bv2R*m3*(m2**2 + m4**2 - tman))*uman**2 - Bv1R*m3*uman**3)*xv) ))/
     &  (sman - mq1c**2))

  
	end


! --- MtBoxV  ---	

	double precision function MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
        
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc
        double precision mq2c,xv
              
          MtBoxV = 0d0
          MtBoxV = dreal((0.0031662869888230555*(A2Lc*CLc*(4*Bs0R*m1*m2**2 + 3*Bs2L*m1**2*m2**2 + 2*Bs3L*m1**2*m2**2 + 4*Bv0L*m1**2*m2**2 - 3*Bv1R*m1**3*m2**2 - 2*Bv1R*m1*m2**4 + 3*Bv2R*m1*m2**4 + 
     &     2*Bv3R*m1*m2**4 - 6*Bs0L*m1**2*m3 - 6*Bs1R*m1**3*m3 + 3*Bv1L*m1**4*m3 - 4*Bs0L*m2**2*m3 - 4*Bs1R*m1*m2**2*m3 + 6*Bs2R*m1*m2**2*m3 + 
     &     4*Bs3R*m1*m2**2*m3 + 6*Bv0R*m1*m2**2*m3 + 5*Bv1L*m1**2*m2**2*m3 - 3*Bv2L*m1**2*m2**2*m3 - 2*Bv3L*m1**2*m2**2*m3 + 2*Bv1L*m2**4*m3 - 
     &     3*Bv2L*m2**4*m3 - 2*Bv3L*m2**4*m3 + 6*Bs0R*m1*m3**2 + 3*Bs3L*m1**2*m3**2 - 3*Bv1R*m1**3*m3**2 + 3*Bs2L*m2**2*m3**2 + 2*Bs3L*m2**2*m3**2 + 
     &     2*Bv0L*m2**2*m3**2 - 2*Bv1R*m1*m2**2*m3**2 + 3*Bv2R*m1*m2**2*m3**2 + 5*Bv3R*m1*m2**2*m3**2 + 6*Bs3R*m1*m3**3 - 3*Bv3L*m1**2*m3**3 - 
     &     3*Bv3L*m2**2*m3**3 + 3*Bs3L*m3**4 + 3*Bv3R*m1*m3**4 + 2*Bs0R*m1*m4**2 - Bs2L*m1**2*m4**2 + Bs3L*m1**2*m4**2 - 2*Bv0L*m1**2*m4**2 - 
     &     Bv1R*m1*m2**2*m4**2 - Bv2R*m1*m2**2*m4**2 + Bv3R*m1*m2**2*m4**2 - 2*Bs0L*m3*m4**2 - 2*Bs1R*m1*m3*m4**2 - 2*Bs2R*m1*m3*m4**2 + 2*Bs3R*m1*m3*m4**2 - 
     &     2*Bv0R*m1*m3*m4**2 + Bv1L*m1**2*m3*m4**2 + Bv2L*m1**2*m3*m4**2 - Bv3L*m1**2*m3*m4**2 + Bv1L*m2**2*m3*m4**2 + Bv2L*m2**2*m3*m4**2 - 
     &     Bv3L*m2**2*m3*m4**2 - Bs2L*m3**2*m4**2 + Bs3L*m3**2*m4**2 - Bv1R*m1*m3**2*m4**2 - Bv2R*m1*m3**2*m4**2 + Bv3R*m1*m3**2*m4**2 - 2*Bs0R*m1*sman - 
     &     Bs3L*m1**2*sman + Bv0L*m1**2*sman + Bv0L*m2**2*sman + 2*Bv1R*m1*m2**2*sman - Bv3R*m1*m2**2*sman + 4*Bs0L*m3*sman + 4*Bs1R*m1*m3*sman - 
     &     2*Bs3R*m1*m3*sman - 5*Bv1L*m1**2*m3*sman + Bv3L*m1**2*m3*sman - 4*Bv1L*m2**2*m3*sman + 3*Bv2L*m2**2*m3*sman + 3*Bv3L*m2**2*m3*sman - 
     &     Bs3L*m3**2*sman + Bv0L*m3**2*sman + 2*Bv1R*m1*m3**2*sman - Bv3R*m1*m3**2*sman + 3*Bv3L*m3**3*sman + Bv0L*m4**2*sman - Bv1L*m3*m4**2*sman - 
     &     Bv2L*m3*m4**2*sman + Bv3L*m3*m4**2*sman - Bv0L*sman**2 + 2*Bv1L*m3*sman**2 - Bv3L*m3*sman**2 + Bs2L*m1**2*tman + Bv0L*m1**2*tman - 
     &     3*Bs2L*m2**2*tman - 2*Bs3L*m2**2*tman - 3*Bv0L*m2**2*tman + Bv2R*m1*m2**2*tman + 2*Bs2R*m1*m3*tman + 2*Bv0R*m1*m3*tman - Bv2L*m1**2*m3*tman - 
     &     Bv2L*m2**2*m3*tman + Bs2L*m3**2*tman - 3*Bs3L*m3**2*tman + Bv0L*m3**2*tman + Bv2R*m1*m3**2*tman + Bs2L*m4**2*tman - Bs3L*m4**2*tman + 
     &     Bv0L*m4**2*tman + Bs3L*sman*tman + Bv2L*m3*sman*tman - Bs2L*tman**2 - Bv0L*tman**2 - 4*Bs0R*m1*uman - 2*Bs3L*m1**2*uman - Bv0L*m1**2*uman + 
     &     3*Bv1R*m1**3*uman - Bv0L*m2**2*uman + 3*Bv1R*m1*m2**2*uman - 3*Bv2R*m1*m2**2*uman - 4*Bv3R*m1*m2**2*uman + 2*Bs0L*m3*uman + 2*Bs1R*m1*m3*uman - 
     &     4*Bs3R*m1*m3*uman - Bv1L*m1**2*m3*uman + 2*Bv3L*m1**2*m3*uman - Bv1L*m2**2*m3*uman + 2*Bv3L*m2**2*m3*uman - 2*Bs3L*m3**2*uman - Bv0L*m3**2*uman + 
     &     Bv1R*m1*m3**2*uman - 5*Bv3R*m1*m3**2*uman - Bv0L*m4**2*uman + Bv1R*m1*m4**2*uman + Bv2R*m1*m4**2*uman - Bv3R*m1*m4**2*uman - 2*Bv1R*m1*sman*uman + 
     &     Bv3R*m1*sman*uman + Bv1L*m3*sman*uman - 2*Bv3L*m3*sman*uman + 2*Bs3L*tman*uman - Bv2R*m1*tman*uman + Bv0L*uman**2 - Bv1R*m1*uman**2 + 
     &     2*Bv3R*m1*uman**2 + (-m2**2 + tman)*(Bv2R*m1*m2**4 + 2*Bs0L*m1**2*m3 + 2*Bs1R*m1**3*m3 - Bv1L*m1**4*m3 + 2*Bs2R*m1*m2**2*m3 + 2*Bv0R*m1*m2**2*m3 - 
     &        Bv1L*m1**2*m2**2*m3 - Bv2L*m1**2*m2**2*m3 - Bv2L*m2**4*m3 - 2*Bs0R*m1*m3**2 - Bs3L*m1**2*m3**2 + 2*Bv0L*m2**2*m3**2 + Bv2R*m1*m2**2*m3**2 - 
     &        Bv3R*m1*m2**2*m3**2 - 2*Bs3R*m1*m3**3 + Bv3L*m1**2*m3**3 + Bv3L*m2**2*m3**3 - Bs3L*m3**4 - Bv3R*m1*m3**4 - 2*Bs0R*m1*m4**2 - Bs3L*m1**2*m4**2 + 
     &        2*Bv0L*m1**2*m4**2 + Bv2R*m1*m2**2*m4**2 - Bv3R*m1*m2**2*m4**2 + 2*Bs0L*m3*m4**2 + 2*Bs1R*m1*m3*m4**2 + 2*Bs2R*m1*m3*m4**2 - 
     &        2*Bs3R*m1*m3*m4**2 + 2*Bv0R*m1*m3*m4**2 - Bv1L*m1**2*m3*m4**2 - Bv2L*m1**2*m3*m4**2 + Bv3L*m1**2*m3*m4**2 - Bv1L*m2**2*m3*m4**2 - 
     &        Bv2L*m2**2*m3*m4**2 + Bv3L*m2**2*m3*m4**2 - Bs3L*m3**2*m4**2 + Bv2R*m1*m3**2*m4**2 - Bv3R*m1*m3**2*m4**2 + 2*Bs0R*m1*sman + Bs3L*m1**2*sman - 
     &        Bv0L*m1**2*sman - Bv0L*m2**2*sman + Bv3R*m1*m2**2*sman + 2*Bs3R*m1*m3*sman + Bv1L*m1**2*m3*sman - Bv3L*m1**2*m3*sman + Bv2L*m2**2*m3*sman - 
     &        Bv3L*m2**2*m3*sman + Bs3L*m3**2*sman - Bv0L*m3**2*sman + Bv3R*m1*m3**2*sman - Bv3L*m3**3*sman - Bv0L*m4**2*sman + Bv1L*m3*m4**2*sman + 
     &        Bv2L*m3*m4**2*sman - Bv3L*m3*m4**2*sman + Bv0L*sman**2 + Bv3L*m3*sman**2 + Bs2L*(m1**2 + m3**2 - tman)*(m2**2 + m4**2 - tman) - 
     &        Bv0L*m1**2*tman - Bv0L*m2**2*tman - Bv2R*m1*m2**2*tman - 2*Bs2R*m1*m3*tman - 2*Bv0R*m1*m3*tman + Bv2L*m1**2*m3*tman + Bv2L*m2**2*m3*tman + 
     &        Bs3L*m3**2*tman - Bv0L*m3**2*tman - Bv2R*m1*m3**2*tman + Bs3L*m4**2*tman - Bv0L*m4**2*tman - Bs3L*sman*tman - Bv2L*m3*sman*tman + Bv0L*tman**2 + 
     &        Bv1R*m1*(m2**2 + m3**2 - uman)*(m1**2 + m4**2 - uman) + 
     &        (-2*Bs0L*m3 - 2*Bs1R*m1*m3 + Bv1L*m1**2*m3 + Bv1L*m2**2*m3 + Bv3R*m1*m3**2 + Bv3R*m1*m4**2 + Bv0L*(m1**2 + m2**2 + m3**2 + m4**2) - 
     &           (Bv3R*m1 + Bv1L*m3)*sman - Bv2R*m1*(m2**2 + m4**2 - tman))*uman - Bv0L*uman**2)*xv + 
     &     Bs1L*(m1**2 + m3**2 - tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sman + uman - (m2**2 - tman)*(m1**2 + m4**2 - uman)*xv)) + 
     &  A2Rc*CLc*(4*Bs0L*m1*m2**2 + 3*Bs2R*m1**2*m2**2 + 2*Bs3R*m1**2*m2**2 + 4*Bv0R*m1**2*m2**2 - 3*Bv1L*m1**3*m2**2 - 2*Bv1L*m1*m2**4 + 3*Bv2L*m1*m2**4 + 
     &     2*Bv3L*m1*m2**4 - 6*Bs0R*m1**2*m3 - 6*Bs1L*m1**3*m3 + 3*Bv1R*m1**4*m3 - 4*Bs0R*m2**2*m3 - 4*Bs1L*m1*m2**2*m3 + 6*Bs2L*m1*m2**2*m3 + 
     &     4*Bs3L*m1*m2**2*m3 + 6*Bv0L*m1*m2**2*m3 + 5*Bv1R*m1**2*m2**2*m3 - 3*Bv2R*m1**2*m2**2*m3 - 2*Bv3R*m1**2*m2**2*m3 + 2*Bv1R*m2**4*m3 - 
     &     3*Bv2R*m2**4*m3 - 2*Bv3R*m2**4*m3 + 6*Bs0L*m1*m3**2 + 3*Bs3R*m1**2*m3**2 - 3*Bv1L*m1**3*m3**2 + 3*Bs2R*m2**2*m3**2 + 2*Bs3R*m2**2*m3**2 + 
     &     2*Bv0R*m2**2*m3**2 - 2*Bv1L*m1*m2**2*m3**2 + 3*Bv2L*m1*m2**2*m3**2 + 5*Bv3L*m1*m2**2*m3**2 + 6*Bs3L*m1*m3**3 - 3*Bv3R*m1**2*m3**3 - 
     &     3*Bv3R*m2**2*m3**3 + 3*Bs3R*m3**4 + 3*Bv3L*m1*m3**4 + 2*Bs0L*m1*m4**2 - Bs2R*m1**2*m4**2 + Bs3R*m1**2*m4**2 - 2*Bv0R*m1**2*m4**2 - 
     &     Bv1L*m1*m2**2*m4**2 - Bv2L*m1*m2**2*m4**2 + Bv3L*m1*m2**2*m4**2 - 2*Bs0R*m3*m4**2 - 2*Bs1L*m1*m3*m4**2 - 2*Bs2L*m1*m3*m4**2 + 2*Bs3L*m1*m3*m4**2 - 
     &     2*Bv0L*m1*m3*m4**2 + Bv1R*m1**2*m3*m4**2 + Bv2R*m1**2*m3*m4**2 - Bv3R*m1**2*m3*m4**2 + Bv1R*m2**2*m3*m4**2 + Bv2R*m2**2*m3*m4**2 - 
     &     Bv3R*m2**2*m3*m4**2 - Bs2R*m3**2*m4**2 + Bs3R*m3**2*m4**2 - Bv1L*m1*m3**2*m4**2 - Bv2L*m1*m3**2*m4**2 + Bv3L*m1*m3**2*m4**2 - 2*Bs0L*m1*sman - 
     &     Bs3R*m1**2*sman + Bv0R*m1**2*sman + Bv0R*m2**2*sman + 2*Bv1L*m1*m2**2*sman - Bv3L*m1*m2**2*sman + 4*Bs0R*m3*sman + 4*Bs1L*m1*m3*sman - 
     &     2*Bs3L*m1*m3*sman - 5*Bv1R*m1**2*m3*sman + Bv3R*m1**2*m3*sman - 4*Bv1R*m2**2*m3*sman + 3*Bv2R*m2**2*m3*sman + 3*Bv3R*m2**2*m3*sman - 
     &     Bs3R*m3**2*sman + Bv0R*m3**2*sman + 2*Bv1L*m1*m3**2*sman - Bv3L*m1*m3**2*sman + 3*Bv3R*m3**3*sman + Bv0R*m4**2*sman - Bv1R*m3*m4**2*sman - 
     &     Bv2R*m3*m4**2*sman + Bv3R*m3*m4**2*sman - Bv0R*sman**2 + 2*Bv1R*m3*sman**2 - Bv3R*m3*sman**2 + Bs2R*m1**2*tman + Bv0R*m1**2*tman - 
     &     3*Bs2R*m2**2*tman - 2*Bs3R*m2**2*tman - 3*Bv0R*m2**2*tman + Bv2L*m1*m2**2*tman + 2*Bs2L*m1*m3*tman + 2*Bv0L*m1*m3*tman - Bv2R*m1**2*m3*tman - 
     &     Bv2R*m2**2*m3*tman + Bs2R*m3**2*tman - 3*Bs3R*m3**2*tman + Bv0R*m3**2*tman + Bv2L*m1*m3**2*tman + Bs2R*m4**2*tman - Bs3R*m4**2*tman + 
     &     Bv0R*m4**2*tman + Bs3R*sman*tman + Bv2R*m3*sman*tman - Bs2R*tman**2 - Bv0R*tman**2 - 4*Bs0L*m1*uman - 2*Bs3R*m1**2*uman - Bv0R*m1**2*uman + 
     &     3*Bv1L*m1**3*uman - Bv0R*m2**2*uman + 3*Bv1L*m1*m2**2*uman - 3*Bv2L*m1*m2**2*uman - 4*Bv3L*m1*m2**2*uman + 2*Bs0R*m3*uman + 2*Bs1L*m1*m3*uman - 
     &     4*Bs3L*m1*m3*uman - Bv1R*m1**2*m3*uman + 2*Bv3R*m1**2*m3*uman - Bv1R*m2**2*m3*uman + 2*Bv3R*m2**2*m3*uman - 2*Bs3R*m3**2*uman - Bv0R*m3**2*uman + 
     &     Bv1L*m1*m3**2*uman - 5*Bv3L*m1*m3**2*uman - Bv0R*m4**2*uman + Bv1L*m1*m4**2*uman + Bv2L*m1*m4**2*uman - Bv3L*m1*m4**2*uman - 2*Bv1L*m1*sman*uman + 
     &     Bv3L*m1*sman*uman + Bv1R*m3*sman*uman - 2*Bv3R*m3*sman*uman + 2*Bs3R*tman*uman - Bv2L*m1*tman*uman + Bv0R*uman**2 - Bv1L*m1*uman**2 + 
     &     2*Bv3L*m1*uman**2 + (-m2**2 + tman)*(Bv2L*m1*m2**4 + 2*Bs0R*m1**2*m3 + 2*Bs1L*m1**3*m3 - Bv1R*m1**4*m3 + 2*Bs2L*m1*m2**2*m3 + 2*Bv0L*m1*m2**2*m3 - 
     &        Bv1R*m1**2*m2**2*m3 - Bv2R*m1**2*m2**2*m3 - Bv2R*m2**4*m3 - 2*Bs0L*m1*m3**2 - Bs3R*m1**2*m3**2 + 2*Bv0R*m2**2*m3**2 + Bv2L*m1*m2**2*m3**2 - 
     &        Bv3L*m1*m2**2*m3**2 - 2*Bs3L*m1*m3**3 + Bv3R*m1**2*m3**3 + Bv3R*m2**2*m3**3 - Bs3R*m3**4 - Bv3L*m1*m3**4 - 2*Bs0L*m1*m4**2 - Bs3R*m1**2*m4**2 + 
     &        2*Bv0R*m1**2*m4**2 + Bv2L*m1*m2**2*m4**2 - Bv3L*m1*m2**2*m4**2 + 2*Bs0R*m3*m4**2 + 2*Bs1L*m1*m3*m4**2 + 2*Bs2L*m1*m3*m4**2 - 
     &        2*Bs3L*m1*m3*m4**2 + 2*Bv0L*m1*m3*m4**2 - Bv1R*m1**2*m3*m4**2 - Bv2R*m1**2*m3*m4**2 + Bv3R*m1**2*m3*m4**2 - Bv1R*m2**2*m3*m4**2 - 
     &        Bv2R*m2**2*m3*m4**2 + Bv3R*m2**2*m3*m4**2 - Bs3R*m3**2*m4**2 + Bv2L*m1*m3**2*m4**2 - Bv3L*m1*m3**2*m4**2 + 2*Bs0L*m1*sman + Bs3R*m1**2*sman - 
     &        Bv0R*m1**2*sman - Bv0R*m2**2*sman + Bv3L*m1*m2**2*sman + 2*Bs3L*m1*m3*sman + Bv1R*m1**2*m3*sman - Bv3R*m1**2*m3*sman + Bv2R*m2**2*m3*sman - 
     &        Bv3R*m2**2*m3*sman + Bs3R*m3**2*sman - Bv0R*m3**2*sman + Bv3L*m1*m3**2*sman - Bv3R*m3**3*sman - Bv0R*m4**2*sman + Bv1R*m3*m4**2*sman + 
     &        Bv2R*m3*m4**2*sman - Bv3R*m3*m4**2*sman + Bv0R*sman**2 + Bv3R*m3*sman**2 + Bs2R*(m1**2 + m3**2 - tman)*(m2**2 + m4**2 - tman) - 
     &        Bv0R*m1**2*tman - Bv0R*m2**2*tman - Bv2L*m1*m2**2*tman - 2*Bs2L*m1*m3*tman - 2*Bv0L*m1*m3*tman + Bv2R*m1**2*m3*tman + Bv2R*m2**2*m3*tman + 
     &        Bs3R*m3**2*tman - Bv0R*m3**2*tman - Bv2L*m1*m3**2*tman + Bs3R*m4**2*tman - Bv0R*m4**2*tman - Bs3R*sman*tman - Bv2R*m3*sman*tman + Bv0R*tman**2 + 
     &        Bv1L*m1*(m2**2 + m3**2 - uman)*(m1**2 + m4**2 - uman) + 
     &        (-2*Bs0R*m3 - 2*Bs1L*m1*m3 + Bv1R*m1**2*m3 + Bv1R*m2**2*m3 + Bv3L*m1*m3**2 + Bv3L*m1*m4**2 + Bv0R*(m1**2 + m2**2 + m3**2 + m4**2) - 
     &           (Bv3L*m1 + Bv1R*m3)*sman - Bv2L*m1*(m2**2 + m4**2 - tman))*uman - Bv0R*uman**2)*xv + 
     &     Bs1R*(m1**2 + m3**2 - tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sman + uman - (m2**2 - tman)*(m1**2 + m4**2 - uman)*xv)) ))/(tman - mq2c**2))
                  

	end


! --- MuBoxV  ---	

	double precision function MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)
        
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc
        double precision mq3c,xv
                        
          MuBoxV = 0d0
          MuBoxV = dreal((0.0015831434944115277*(A3Rc*DLc*mq3c*(-8*Bs0L*m1**2 + 4*Bs1R*m1**3 + 4*Bs2R*m1*m2**2 - 8*Bv0R*m1*m2**2 + 4*Bv1L*m1**2*m2**2 - 8*Bv2L*m1**2*m2**2 - 4*Bv2L*m2**4 + 
     &     32*Bs0R*m1*m3 + 8*Bs1L*m1**2*m3 - 4*Bs2L*m1**2*m3 + 4*Bs3L*m1**2*m3 - 16*Bv0L*m1**2*m3 - 4*Bv1R*m1**3*m3 - 4*Bs2L*m2**2*m3 - 16*Bv0L*m2**2*m3 - 
     &     4*Bv1R*m1*m2**2*m3 + 8*Bv2R*m1*m2**2*m3 + 4*Bv3R*m1*m2**2*m3 - 8*Bs0L*m3**2 + 4*Bs1R*m1*m3**2 + 4*Bs2R*m1*m3**2 + 8*Bs3R*m1*m3**2 - 
     &     8*Bv0R*m1*m3**2 + 4*Bv1L*m1**2*m3**2 - 4*Bv2L*m1**2*m3**2 - 4*Bv3L*m1**2*m3**2 - 8*Bv2L*m2**2*m3**2 - 4*Bv3L*m2**2*m3**2 + 4*Bs3L*m3**3 + 
     &     4*Bv3R*m1*m3**3 + 4*Bv2L*m2**2*sman + 4*Bs2L*m3*sman + 16*Bv0L*m3*sman + 4*Bv1R*m1*m3*sman + 4*Bv2L*m3**2*sman + 4*Bv3L*m3**2*sman + 8*Bs0L*tman - 
     &     4*Bs1R*m1*tman + 4*Bv2L*m2**2*tman - 4*Bs3L*m3*tman - 4*Bs2R*m1*uman + 8*Bv0R*m1*uman - 4*Bv1L*m1**2*uman + 4*Bv2L*m1**2*uman + 4*Bv2L*m2**2*uman - 
     &     4*Bv3R*m1*m3*uman - 4*Bv2L*sman*uman - (2*Bv1R*m1**3*m2**2*m3 + 2*Bv2R*m1*m2**4*m3 - 4*Bs0L*m1**2*m3**2 - 2*Bs1R*m1**3*m3**2 + 2*Bv1L*m1**4*m3**2 - 
     &        2*Bs2R*m1*m2**2*m3**2 - 4*Bv0R*m1*m2**2*m3**2 + 2*Bv2L*m1**2*m2**2*m3**2 - 2*Bs3L*m1**2*m3**3 - 2*Bv3R*m1*m2**2*m3**3 + 2*Bs3R*m1*m3**4 - 
     &        2*Bv3L*m1**2*m3**4 - 8*Bs0L*m1**2*m4**2 - 2*Bs1R*m1**3*m4**2 - 2*Bs2R*m1*m2**2*m4**2 - 8*Bv0R*m1*m2**2*m4**2 + 2*Bv1L*m1**2*m2**2*m4**2 + 
     &        2*Bv2L*m2**4*m4**2 + 8*Bs0R*m1*m3*m4**2 - 2*Bs3L*m1**2*m3*m4**2 - 4*Bv0L*m1**2*m3*m4**2 + 2*Bv1R*m1**3*m3*m4**2 - 4*Bv0L*m2**2*m3*m4**2 + 
     &        2*Bv1R*m1*m2**2*m3*m4**2 + 4*Bv2R*m1*m2**2*m3*m4**2 - 2*Bv3R*m1*m2**2*m3*m4**2 - 8*Bs0L*m3**2*m4**2 - 2*Bs1R*m1*m3**2*m4**2 - 
     &        2*Bs2R*m1*m3**2*m4**2 + 4*Bs3R*m1*m3**2*m4**2 - 8*Bv0R*m1*m3**2*m4**2 + 2*Bv1L*m1**2*m3**2*m4**2 + 2*Bv2L*m1**2*m3**2*m4**2 - 
     &        2*Bv3L*m1**2*m3**2*m4**2 - 2*Bv3L*m2**2*m3**2*m4**2 - 2*Bs3L*m3**3*m4**2 - 2*Bv3R*m1*m3**3*m4**2 - 4*Bs0L*m4**4 - 2*Bs1R*m1*m4**4 - 
     &        2*Bs2R*m1*m4**4 + 2*Bs3R*m1*m4**4 - 4*Bv0R*m1*m4**4 + 2*Bv1L*m2**2*m4**4 + 2*Bv2L*m2**2*m4**4 - 2*Bv3L*m2**2*m4**4 - 2*Bs3L*m3*m4**4 + 
     &        2*Bv1R*m1*m3*m4**4 + 2*Bv2R*m1*m3*m4**4 - 2*Bv3R*m1*m3*m4**4 + 4*Bs0L*m1**2*sman + 2*Bs1R*m1**3*sman - Bv1L*m1**4*sman + 2*Bs2R*m1*m2**2*sman + 
     &        4*Bv0R*m1*m2**2*sman - Bv1L*m1**2*m2**2*sman - Bv2L*m1**2*m2**2*sman - Bv2L*m2**4*sman + 2*Bs3L*m1**2*m3*sman + 2*Bv3R*m1*m2**2*m3*sman - 
     &        4*Bs3R*m1*m3**2*sman - Bv1L*m1**2*m3**2*sman + 3*Bv3L*m1**2*m3**2*sman - Bv2L*m2**2*m3**2*sman + Bv3L*m2**2*m3**2*sman + Bv3L*m3**4*sman + 
     &        4*Bs0L*m4**2*sman + 2*Bs1R*m1*m4**2*sman + 2*Bs2R*m1*m4**2*sman - 4*Bs3R*m1*m4**2*sman + 4*Bv0R*m1*m4**2*sman - 2*Bv1L*m1**2*m4**2*sman - 
     &        Bv2L*m1**2*m4**2*sman + Bv3L*m1**2*m4**2*sman - Bv1L*m2**2*m4**2*sman - 2*Bv2L*m2**2*m4**2*sman + 3*Bv3L*m2**2*m4**2*sman + 
     &        2*Bs3L*m3*m4**2*sman + 4*Bv0L*m3*m4**2*sman + 2*Bv3R*m1*m3*m4**2*sman - Bv1L*m3**2*m4**2*sman - Bv2L*m3**2*m4**2*sman + 
     &        2*Bv3L*m3**2*m4**2*sman - Bv1L*m4**4*sman - Bv2L*m4**4*sman + Bv3L*m4**4*sman + 2*Bs3R*m1*sman**2 + Bv1L*m1**2*sman**2 - Bv3L*m1**2*sman**2 + 
     &        Bv2L*m2**2*sman**2 - Bv3L*m2**2*sman**2 - 2*Bv3L*m3**2*sman**2 + Bv1L*m4**2*sman**2 + Bv2L*m4**2*sman**2 - 2*Bv3L*m4**2*sman**2 + Bv3L*sman**3 + 
     &        Bv1L*m1**4*tman + Bv1L*m1**2*m2**2*tman + Bv2L*m1**2*m2**2*tman + Bv2L*m2**4*tman - 2*Bv1R*m1**3*m3*tman - 4*Bv2R*m1*m2**2*m3*tman + 
     &        2*Bs2R*m1*m3**2*tman + 4*Bv0R*m1*m3**2*tman + Bv1L*m1**2*m3**2*tman - 2*Bv2L*m1**2*m3**2*tman - Bv3L*m1**2*m3**2*tman + Bv2L*m2**2*m3**2*tman - 
     &        Bv3L*m2**2*m3**2*tman + 2*Bv3R*m1*m3**3*tman - Bv3L*m3**4*tman + 4*Bs0L*m4**2*tman + 2*Bs2R*m1*m4**2*tman + 4*Bv0R*m1*m4**2*tman + 
     &        2*Bv1L*m1**2*m4**2*tman + Bv2L*m1**2*m4**2*tman - Bv3L*m1**2*m4**2*tman + Bv1L*m2**2*m4**2*tman - Bv3L*m2**2*m4**2*tman - 
     &        2*Bv1R*m1*m3*m4**2*tman - 4*Bv2R*m1*m3*m4**2*tman + 2*Bv3R*m1*m3*m4**2*tman + Bv1L*m3**2*m4**2*tman + Bv2L*m3**2*m4**2*tman - 
     &        2*Bv3L*m3**2*m4**2*tman + Bv1L*m4**4*tman + Bv2L*m4**4*tman - Bv3L*m4**4*tman - 2*Bs2R*m1*sman*tman - 4*Bv0R*m1*sman*tman + 
     &        Bv2L*m1**2*sman*tman + Bv3L*m1**2*sman*tman + Bv2L*m2**2*sman*tman + Bv3L*m2**2*sman*tman - 2*Bv3R*m1*m3*sman*tman + Bv2L*m3**2*sman*tman + 
     &        Bv3L*m3**2*sman*tman + Bv2L*m4**2*sman*tman + Bv3L*m4**2*sman*tman - Bv2L*sman**2*tman - Bv1L*m1**2*tman**2 - Bv2L*m1**2*tman**2 - 
     &        2*Bv2L*m2**2*tman**2 + 2*Bv2R*m1*m3*tman**2 - Bv2L*m3**2*tman**2 + Bv3L*m3**2*tman**2 - Bv1L*m4**2*tman**2 - 2*Bv2L*m4**2*tman**2 + 
     &        Bv3L*m4**2*tman**2 - Bv3L*sman*tman**2 + Bv2L*tman**3 + 2*Bs2L*m3*(m2**2 + m4**2 - tman)*(m1**2 + m4**2 - uman) + 
     &        2*Bs1L*m3*(m1**2 + m4**2 - uman)**2 - (2*Bv1R*m1*m2**2*m3 - 4*Bs0L*m3**2 - 2*Bs1R*m1*m3**2 - Bv3L*m1**2*m3**2 - Bv3L*m2**2*m3**2 - 
     &           2*Bs3L*m3**3 - Bv3L*m3**4 - 4*Bs0L*m4**2 - 2*Bs1R*m1*m4**2 - 4*Bv0R*m1*m4**2 - Bv3L*m1**2*m4**2 - Bv3L*m2**2*m4**2 - 2*Bs3L*m3*m4**2 + 
     &           2*Bv1R*m1*m3*m4**2 - 2*Bv3L*m3**2*m4**2 - Bv3L*m4**4 + 4*Bs0L*sman + 2*Bs1R*m1*sman + Bv3L*m1**2*sman + Bv3L*m2**2*sman + 2*Bs3L*m3*sman + 
     &           Bv3L*m3**2*sman + Bv3L*m4**2*sman + Bv2L*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1R*m1*m3*tman + 
     &           Bv1L*(m1**4 + m3**2*m4**2 + m4**4 - m3**2*sman - m4**2*sman + sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + 3*m3**2 + 2*m4**2 - sman + tman) + m2**2*(3*m4**2 - sman + tman)))*uman + 
     &        (Bv1L*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - Bv3L*(m3**2 + m4**2 - sman) + Bv2L*(m2**2 + m4**2 - tman))*uman**2 - Bv1L*uman**3)*xv) + 
     &  A3Lc*DRc*mq3c*(-8*Bs0R*m1**2 + 4*Bs1L*m1**3 + 4*Bs2L*m1*m2**2 - 8*Bv0L*m1*m2**2 + 4*Bv1R*m1**2*m2**2 - 8*Bv2R*m1**2*m2**2 - 4*Bv2R*m2**4 + 
     &     32*Bs0L*m1*m3 + 8*Bs1R*m1**2*m3 - 4*Bs2R*m1**2*m3 + 4*Bs3R*m1**2*m3 - 16*Bv0R*m1**2*m3 - 4*Bv1L*m1**3*m3 - 4*Bs2R*m2**2*m3 - 16*Bv0R*m2**2*m3 - 
     &     4*Bv1L*m1*m2**2*m3 + 8*Bv2L*m1*m2**2*m3 + 4*Bv3L*m1*m2**2*m3 - 8*Bs0R*m3**2 + 4*Bs1L*m1*m3**2 + 4*Bs2L*m1*m3**2 + 8*Bs3L*m1*m3**2 - 
     &     8*Bv0L*m1*m3**2 + 4*Bv1R*m1**2*m3**2 - 4*Bv2R*m1**2*m3**2 - 4*Bv3R*m1**2*m3**2 - 8*Bv2R*m2**2*m3**2 - 4*Bv3R*m2**2*m3**2 + 4*Bs3R*m3**3 + 
     &     4*Bv3L*m1*m3**3 + 4*Bv2R*m2**2*sman + 4*Bs2R*m3*sman + 16*Bv0R*m3*sman + 4*Bv1L*m1*m3*sman + 4*Bv2R*m3**2*sman + 4*Bv3R*m3**2*sman + 8*Bs0R*tman - 
     &     4*Bs1L*m1*tman + 4*Bv2R*m2**2*tman - 4*Bs3R*m3*tman - 4*Bs2L*m1*uman + 8*Bv0L*m1*uman - 4*Bv1R*m1**2*uman + 4*Bv2R*m1**2*uman + 4*Bv2R*m2**2*uman - 
     &     4*Bv3L*m1*m3*uman - 4*Bv2R*sman*uman - (2*Bv1L*m1**3*m2**2*m3 + 2*Bv2L*m1*m2**4*m3 - 4*Bs0R*m1**2*m3**2 - 2*Bs1L*m1**3*m3**2 + 2*Bv1R*m1**4*m3**2 - 
     &        2*Bs2L*m1*m2**2*m3**2 - 4*Bv0L*m1*m2**2*m3**2 + 2*Bv2R*m1**2*m2**2*m3**2 - 2*Bs3R*m1**2*m3**3 - 2*Bv3L*m1*m2**2*m3**3 + 2*Bs3L*m1*m3**4 - 
     &        2*Bv3R*m1**2*m3**4 - 8*Bs0R*m1**2*m4**2 - 2*Bs1L*m1**3*m4**2 - 2*Bs2L*m1*m2**2*m4**2 - 8*Bv0L*m1*m2**2*m4**2 + 2*Bv1R*m1**2*m2**2*m4**2 + 
     &        2*Bv2R*m2**4*m4**2 + 8*Bs0L*m1*m3*m4**2 - 2*Bs3R*m1**2*m3*m4**2 - 4*Bv0R*m1**2*m3*m4**2 + 2*Bv1L*m1**3*m3*m4**2 - 4*Bv0R*m2**2*m3*m4**2 + 
     &        2*Bv1L*m1*m2**2*m3*m4**2 + 4*Bv2L*m1*m2**2*m3*m4**2 - 2*Bv3L*m1*m2**2*m3*m4**2 - 8*Bs0R*m3**2*m4**2 - 2*Bs1L*m1*m3**2*m4**2 - 
     &        2*Bs2L*m1*m3**2*m4**2 + 4*Bs3L*m1*m3**2*m4**2 - 8*Bv0L*m1*m3**2*m4**2 + 2*Bv1R*m1**2*m3**2*m4**2 + 2*Bv2R*m1**2*m3**2*m4**2 - 
     &        2*Bv3R*m1**2*m3**2*m4**2 - 2*Bv3R*m2**2*m3**2*m4**2 - 2*Bs3R*m3**3*m4**2 - 2*Bv3L*m1*m3**3*m4**2 - 4*Bs0R*m4**4 - 2*Bs1L*m1*m4**4 - 
     &        2*Bs2L*m1*m4**4 + 2*Bs3L*m1*m4**4 - 4*Bv0L*m1*m4**4 + 2*Bv1R*m2**2*m4**4 + 2*Bv2R*m2**2*m4**4 - 2*Bv3R*m2**2*m4**4 - 2*Bs3R*m3*m4**4 + 
     &        2*Bv1L*m1*m3*m4**4 + 2*Bv2L*m1*m3*m4**4 - 2*Bv3L*m1*m3*m4**4 + 4*Bs0R*m1**2*sman + 2*Bs1L*m1**3*sman - Bv1R*m1**4*sman + 2*Bs2L*m1*m2**2*sman + 
     &        4*Bv0L*m1*m2**2*sman - Bv1R*m1**2*m2**2*sman - Bv2R*m1**2*m2**2*sman - Bv2R*m2**4*sman + 2*Bs3R*m1**2*m3*sman + 2*Bv3L*m1*m2**2*m3*sman - 
     &        4*Bs3L*m1*m3**2*sman - Bv1R*m1**2*m3**2*sman + 3*Bv3R*m1**2*m3**2*sman - Bv2R*m2**2*m3**2*sman + Bv3R*m2**2*m3**2*sman + Bv3R*m3**4*sman + 
     &        4*Bs0R*m4**2*sman + 2*Bs1L*m1*m4**2*sman + 2*Bs2L*m1*m4**2*sman - 4*Bs3L*m1*m4**2*sman + 4*Bv0L*m1*m4**2*sman - 2*Bv1R*m1**2*m4**2*sman - 
     &        Bv2R*m1**2*m4**2*sman + Bv3R*m1**2*m4**2*sman - Bv1R*m2**2*m4**2*sman - 2*Bv2R*m2**2*m4**2*sman + 3*Bv3R*m2**2*m4**2*sman + 
     &        2*Bs3R*m3*m4**2*sman + 4*Bv0R*m3*m4**2*sman + 2*Bv3L*m1*m3*m4**2*sman - Bv1R*m3**2*m4**2*sman - Bv2R*m3**2*m4**2*sman + 
     &        2*Bv3R*m3**2*m4**2*sman - Bv1R*m4**4*sman - Bv2R*m4**4*sman + Bv3R*m4**4*sman + 2*Bs3L*m1*sman**2 + Bv1R*m1**2*sman**2 - Bv3R*m1**2*sman**2 + 
     &        Bv2R*m2**2*sman**2 - Bv3R*m2**2*sman**2 - 2*Bv3R*m3**2*sman**2 + Bv1R*m4**2*sman**2 + Bv2R*m4**2*sman**2 - 2*Bv3R*m4**2*sman**2 + Bv3R*sman**3 + 
     &        Bv1R*m1**4*tman + Bv1R*m1**2*m2**2*tman + Bv2R*m1**2*m2**2*tman + Bv2R*m2**4*tman - 2*Bv1L*m1**3*m3*tman - 4*Bv2L*m1*m2**2*m3*tman + 
     &        2*Bs2L*m1*m3**2*tman + 4*Bv0L*m1*m3**2*tman + Bv1R*m1**2*m3**2*tman - 2*Bv2R*m1**2*m3**2*tman - Bv3R*m1**2*m3**2*tman + Bv2R*m2**2*m3**2*tman - 
     &        Bv3R*m2**2*m3**2*tman + 2*Bv3L*m1*m3**3*tman - Bv3R*m3**4*tman + 4*Bs0R*m4**2*tman + 2*Bs2L*m1*m4**2*tman + 4*Bv0L*m1*m4**2*tman + 
     &        2*Bv1R*m1**2*m4**2*tman + Bv2R*m1**2*m4**2*tman - Bv3R*m1**2*m4**2*tman + Bv1R*m2**2*m4**2*tman - Bv3R*m2**2*m4**2*tman - 
     &        2*Bv1L*m1*m3*m4**2*tman - 4*Bv2L*m1*m3*m4**2*tman + 2*Bv3L*m1*m3*m4**2*tman + Bv1R*m3**2*m4**2*tman + Bv2R*m3**2*m4**2*tman - 
     &        2*Bv3R*m3**2*m4**2*tman + Bv1R*m4**4*tman + Bv2R*m4**4*tman - Bv3R*m4**4*tman - 2*Bs2L*m1*sman*tman - 4*Bv0L*m1*sman*tman + 
     &        Bv2R*m1**2*sman*tman + Bv3R*m1**2*sman*tman + Bv2R*m2**2*sman*tman + Bv3R*m2**2*sman*tman - 2*Bv3L*m1*m3*sman*tman + Bv2R*m3**2*sman*tman + 
     &        Bv3R*m3**2*sman*tman + Bv2R*m4**2*sman*tman + Bv3R*m4**2*sman*tman - Bv2R*sman**2*tman - Bv1R*m1**2*tman**2 - Bv2R*m1**2*tman**2 - 
     &        2*Bv2R*m2**2*tman**2 + 2*Bv2L*m1*m3*tman**2 - Bv2R*m3**2*tman**2 + Bv3R*m3**2*tman**2 - Bv1R*m4**2*tman**2 - 2*Bv2R*m4**2*tman**2 + 
     &        Bv3R*m4**2*tman**2 - Bv3R*sman*tman**2 + Bv2R*tman**3 + 2*Bs2R*m3*(m2**2 + m4**2 - tman)*(m1**2 + m4**2 - uman) + 
     &        2*Bs1R*m3*(m1**2 + m4**2 - uman)**2 - (2*Bv1L*m1*m2**2*m3 - 4*Bs0R*m3**2 - 2*Bs1L*m1*m3**2 - Bv3R*m1**2*m3**2 - Bv3R*m2**2*m3**2 - 
     &           2*Bs3R*m3**3 - Bv3R*m3**4 - 4*Bs0R*m4**2 - 2*Bs1L*m1*m4**2 - 4*Bv0L*m1*m4**2 - Bv3R*m1**2*m4**2 - Bv3R*m2**2*m4**2 - 2*Bs3R*m3*m4**2 + 
     &           2*Bv1L*m1*m3*m4**2 - 2*Bv3R*m3**2*m4**2 - Bv3R*m4**4 + 4*Bs0R*sman + 2*Bs1L*m1*sman + Bv3R*m1**2*sman + Bv3R*m2**2*sman + 2*Bs3R*m3*sman + 
     &           Bv3R*m3**2*sman + Bv3R*m4**2*sman + Bv2R*(m1**2 + m2**2 + m3**2 + m4**2)*(m2**2 + m4**2 - tman) - 2*Bv1L*m1*m3*tman + 
     &           Bv1R*(m1**4 + m3**2*m4**2 + m4**4 - m3**2*sman - m4**2*sman + sman**2 + (m3**2 + m4**2)*tman - tman**2 + 
     &              m1**2*(m2**2 + 3*m3**2 + 2*m4**2 - sman + tman) + m2**2*(3*m4**2 - sman + tman)))*uman + 
     &        (Bv1R*(2*m1**2 + m2**2 + m3**2 + 2*m4**2) - Bv3R*(m3**2 + m4**2 - sman) + Bv2R*(m2**2 + m4**2 - tman))*uman**2 - Bv1R*uman**3)*xv) + 
     &  A3Lc*DLc*(-4*Bs1L*m1**2*m2**2 - 4*Bs2L*m1**2*m2**2 - 4*Bs3L*m1**2*m2**2 + 8*Bv0L*m1**2*m2**2 + 8*Bv1R*m1**3*m2**2 + 8*Bv0L*m2**4 + 4*Bv1R*m1*m2**4 - 
     &     4*Bv2R*m1*m2**4 - 4*Bv3R*m1*m2**4 - 16*Bs0L*m1**2*m3 + 8*Bs1R*m1**3*m3 - 8*Bs0L*m2**2*m3 + 4*Bs1R*m1*m2**2*m3 - 4*Bs2R*m1*m2**2*m3 - 
     &     4*Bs3R*m1*m2**2*m3 + 8*Bv0R*m1*m2**2*m3 - 4*Bv1L*m1**2*m2**2*m3 - 4*Bv2L*m1**2*m2**2*m3 - 4*Bv3L*m1**2*m2**2*m3 + 4*Bs1L*m1**2*m3**2 - 
     &     4*Bs2L*m1**2*m3**2 - 4*Bs3L*m1**2*m3**2 - 8*Bv0L*m1**2*m3**2 - 8*Bs2L*m2**2*m3**2 - 8*Bs3L*m2**2*m3**2 - 8*Bv0L*m2**2*m3**2 + 
     &     4*Bv1R*m1*m2**2*m3**2 + 4*Bv2R*m1*m2**2*m3**2 + 4*Bv3R*m1*m2**2*m3**2 - 8*Bs0L*m3**3 + 4*Bs1R*m1*m3**3 + 4*Bs2R*m1*m3**3 + 4*Bs3R*m1*m3**3 - 
     &     8*Bv0R*m1*m3**3 + 4*Bv1L*m1**2*m3**3 - 4*Bv2L*m1**2*m3**3 - 4*Bv3L*m1**2*m3**3 - 8*Bv2L*m2**2*m3**3 - 8*Bv3L*m2**2*m3**3 - 8*Bv0L*m2**2*sman - 
     &     4*Bv1R*m1*m2**2*sman + 8*Bs0L*m3*sman - 4*Bs1R*m1*m3*sman + 4*Bs2L*m3**2*sman + 4*Bs3L*m3**2*sman + 8*Bv0L*m3**2*sman + 4*Bv2L*m3**3*sman + 
     &     4*Bv3L*m3**3*sman + 4*Bs2L*m2**2*tman + 4*Bs3L*m2**2*tman - 4*Bv1R*m1*m2**2*tman + 8*Bs0L*m3*tman - 4*Bs1R*m1*m3*tman + 4*Bv2L*m2**2*m3*tman + 
     &     4*Bv3L*m2**2*m3*tman + 4*Bs1L*m1**2*uman + 4*Bs3L*m1**2*uman - 8*Bv0L*m1**2*uman - 4*Bv1R*m1**3*uman - 8*Bv0L*m2**2*uman - 4*Bv1R*m1*m2**2*uman + 
     &     4*Bv2R*m1*m2**2*uman + 8*Bv3R*m1*m2**2*uman - 4*Bs2R*m1*m3*uman + 4*Bs3R*m1*m3*uman + 8*Bv0R*m1*m3*uman - 4*Bv1L*m1**2*m3*uman + 
     &     4*Bv2L*m1**2*m3*uman + 4*Bv2L*m2**2*m3*uman + 4*Bs3L*m3**2*uman + 4*Bv3R*m1*m3**2*uman + 8*Bv0L*sman*uman + 4*Bv1R*m1*sman*uman - 
     &     4*Bv2L*m3*sman*uman - 4*Bs3L*tman*uman - 4*Bv3R*m1*uman**2 + 16*Bs0R*m1*(-m2**2 + m3**2 + uman) + 
     &     (2*Bv1R*m1**3*m2**4 + 2*Bv2R*m1*m2**6 + 4*Bs0L*m1**2*m2**2*m3 + 2*Bs1R*m1**3*m2**2*m3 + 2*Bv1L*m1**4*m2**2*m3 + 2*Bs2R*m1*m2**4*m3 + 
     &        4*Bv0R*m1*m2**4*m3 + 2*Bv2L*m1**2*m2**4*m3 - 2*Bs3L*m1**2*m2**2*m3**2 - 8*Bv0L*m1**2*m2**2*m3**2 + 2*Bv1R*m1**3*m2**2*m3**2 + 
     &        2*Bv2R*m1*m2**4*m3**2 - 2*Bv3R*m1*m2**4*m3**2 + 4*Bs0L*m1**2*m3**3 + 2*Bs1R*m1**3*m3**3 - 2*Bv1L*m1**4*m3**3 + 2*Bs2R*m1*m2**2*m3**3 - 
     &        2*Bs3R*m1*m2**2*m3**3 + 4*Bv0R*m1*m2**2*m3**3 - 2*Bv2L*m1**2*m2**2*m3**3 - 2*Bv3L*m1**2*m2**2*m3**3 + 2*Bs3L*m1**2*m3**4 - 
     &        2*Bv3R*m1*m2**2*m3**4 - 2*Bs3R*m1*m3**5 + 2*Bv3L*m1**2*m3**5 + 4*Bs0R*m1*m2**2*m4**2 - 2*Bs3L*m1**2*m2**2*m4**2 - 12*Bv0L*m1**2*m2**2*m4**2 + 
     &        4*Bv1R*m1**3*m2**2*m4**2 - 8*Bv0L*m2**4*m4**2 + 2*Bv1R*m1*m2**4*m4**2 + 6*Bv2R*m1*m2**4*m4**2 - 2*Bv3R*m1*m2**4*m4**2 + 16*Bs0L*m1**2*m3*m4**2 + 
     &        4*Bs1R*m1**3*m3*m4**2 + 8*Bs0L*m2**2*m3*m4**2 + 2*Bs1R*m1*m2**2*m3*m4**2 + 6*Bs2R*m1*m2**2*m3*m4**2 - 2*Bs3R*m1*m2**2*m3*m4**2 + 
     &        8*Bv0R*m1*m2**2*m3*m4**2 + 2*Bv1L*m1**2*m2**2*m3*m4**2 + 2*Bv2L*m1**2*m2**2*m3*m4**2 - 2*Bv3L*m1**2*m2**2*m3*m4**2 - 4*Bs0R*m1*m3**2*m4**2 + 
     &        2*Bs3L*m1**2*m3**2*m4**2 - 4*Bv0L*m1**2*m3**2*m4**2 - 8*Bv0L*m2**2*m3**2*m4**2 + 2*Bv1R*m1*m2**2*m3**2*m4**2 + 2*Bv2R*m1*m2**2*m3**2*m4**2 - 
     &        6*Bv3R*m1*m2**2*m3**2*m4**2 + 8*Bs0L*m3**3*m4**2 + 2*Bs1R*m1*m3**3*m4**2 + 2*Bs2R*m1*m3**3*m4**2 - 6*Bs3R*m1*m3**3*m4**2 + 
     &        8*Bv0R*m1*m3**3*m4**2 - 2*Bv1L*m1**2*m3**3*m4**2 - 2*Bv2L*m1**2*m3**3*m4**2 + 2*Bv3L*m1**2*m3**3*m4**2 - 8*Bv0L*m2**2*m4**4 + 
     &        4*Bv1R*m1*m2**2*m4**4 + 4*Bv2R*m1*m2**2*m4**4 - 4*Bv3R*m1*m2**2*m4**4 + 8*Bs0L*m3*m4**4 + 4*Bs1R*m1*m3*m4**4 + 4*Bs2R*m1*m3*m4**4 - 
     &        4*Bs3R*m1*m3*m4**4 + 8*Bv0R*m1*m3*m4**4 + 2*Bs3L*m1**2*m2**2*sman + 6*Bv0L*m1**2*m2**2*sman - 2*Bv1R*m1**3*m2**2*sman + 2*Bv0L*m2**4*sman - 
     &        2*Bv2R*m1*m2**4*sman + 2*Bv3R*m1*m2**4*sman - 4*Bs0L*m1**2*m3*sman - 2*Bs1R*m1**3*m3*sman + Bv1L*m1**4*m3*sman - 2*Bs2R*m1*m2**2*m3*sman + 
     &        2*Bs3R*m1*m2**2*m3*sman - 4*Bv0R*m1*m2**2*m3*sman + Bv1L*m1**2*m2**2*m3*sman + Bv2L*m1**2*m2**2*m3*sman + 2*Bv3L*m1**2*m2**2*m3*sman + 
     &        Bv2L*m2**4*m3*sman - 3*Bs3L*m1**2*m3**2*sman - Bs3L*m2**2*m3**2*sman + 2*Bv0L*m2**2*m3**2*sman + 4*Bv3R*m1*m2**2*m3**2*sman + 
     &        4*Bs3R*m1*m3**3*sman + Bv1L*m1**2*m3**3*sman - 3*Bv3L*m1**2*m3**3*sman + Bv2L*m2**2*m3**3*sman - Bv3L*m2**2*m3**3*sman - Bs3L*m3**4*sman - 
     &        Bv3L*m3**5*sman - Bs3L*m1**2*m4**2*sman + 2*Bv0L*m1**2*m4**2*sman - Bs3L*m2**2*m4**2*sman + 12*Bv0L*m2**2*m4**2*sman - 
     &        2*Bv1R*m1*m2**2*m4**2*sman - 2*Bv2R*m1*m2**2*m4**2*sman + 6*Bv3R*m1*m2**2*m4**2*sman - 8*Bs0L*m3*m4**2*sman - 2*Bs1R*m1*m3*m4**2*sman - 
     &        2*Bs2R*m1*m3*m4**2*sman + 6*Bs3R*m1*m3*m4**2*sman - 4*Bv0R*m1*m3*m4**2*sman + 2*Bv1L*m1**2*m3*m4**2*sman + Bv2L*m1**2*m3*m4**2*sman - 
     &        Bv3L*m1**2*m3*m4**2*sman + Bv1L*m2**2*m3*m4**2*sman + 2*Bv2L*m2**2*m3*m4**2*sman - Bv3L*m2**2*m3*m4**2*sman - 2*Bs3L*m3**2*m4**2*sman + 
     &        2*Bv0L*m3**2*m4**2*sman + Bv1L*m3**3*m4**2*sman + Bv2L*m3**3*m4**2*sman - 2*Bv3L*m3**3*m4**2*sman - Bs3L*m4**4*sman + 2*Bv0L*m4**4*sman + 
     &        Bv1L*m3*m4**4*sman + Bv2L*m3*m4**4*sman - Bv3L*m3*m4**4*sman + Bs3L*m1**2*sman**2 + Bs3L*m2**2*sman**2 - 2*Bv0L*m2**2*sman**2 - 
     &        2*Bv3R*m1*m2**2*sman**2 - 2*Bs3R*m1*m3*sman**2 - Bv1L*m1**2*m3*sman**2 + Bv3L*m1**2*m3*sman**2 - Bv2L*m2**2*m3*sman**2 + Bv3L*m2**2*m3*sman**2 + 
     &        2*Bs3L*m3**2*sman**2 + 2*Bv3L*m3**3*sman**2 + 2*Bs3L*m4**2*sman**2 - 2*Bv0L*m4**2*sman**2 - Bv1L*m3*m4**2*sman**2 - Bv2L*m3*m4**2*sman**2 + 
     &        2*Bv3L*m3*m4**2*sman**2 - Bs3L*sman**3 - Bv3L*m3*sman**3 - 2*Bv0L*m1**2*m2**2*tman - 2*Bv1R*m1**3*m2**2*tman - 2*Bv0L*m2**4*tman - 
     &        4*Bv2R*m1*m2**4*tman - 4*Bs0L*m1**2*m3*tman - 2*Bs1R*m1**3*m3*tman - Bv1L*m1**4*m3*tman - 4*Bs2R*m1*m2**2*m3*tman - 8*Bv0R*m1*m2**2*m3*tman - 
     &        Bv1L*m1**2*m2**2*m3*tman - 3*Bv2L*m1**2*m2**2*m3*tman - Bv2L*m2**4*m3*tman + Bs3L*m1**2*m3**2*tman + 4*Bv0L*m1**2*m3**2*tman + 
     &        Bs3L*m2**2*m3**2*tman - 2*Bv0L*m2**2*m3**2*tman - 2*Bv2R*m1*m2**2*m3**2*tman + 2*Bv3R*m1*m2**2*m3**2*tman - 2*Bs2R*m1*m3**3*tman + 
     &        2*Bs3R*m1*m3**3*tman - 4*Bv0R*m1*m3**3*tman - Bv1L*m1**2*m3**3*tman + 2*Bv2L*m1**2*m3**3*tman + Bv3L*m1**2*m3**3*tman - Bv2L*m2**2*m3**3*tman + 
     &        Bv3L*m2**2*m3**3*tman + Bs3L*m3**4*tman + Bv3L*m3**5*tman + Bs3L*m1**2*m4**2*tman - 2*Bv0L*m1**2*m4**2*tman + Bs3L*m2**2*m4**2*tman + 
     &        4*Bv0L*m2**2*m4**2*tman - 2*Bv1R*m1*m2**2*m4**2*tman - 6*Bv2R*m1*m2**2*m4**2*tman + 2*Bv3R*m1*m2**2*m4**2*tman - 8*Bs0L*m3*m4**2*tman - 
     &        2*Bs1R*m1*m3*m4**2*tman - 6*Bs2R*m1*m3*m4**2*tman + 2*Bs3R*m1*m3*m4**2*tman - 12*Bv0R*m1*m3*m4**2*tman - 2*Bv1L*m1**2*m3*m4**2*tman - 
     &        Bv2L*m1**2*m3*m4**2*tman + Bv3L*m1**2*m3*m4**2*tman - Bv1L*m2**2*m3*m4**2*tman - 2*Bv2L*m2**2*m3*m4**2*tman + Bv3L*m2**2*m3*m4**2*tman + 
     &        2*Bs3L*m3**2*m4**2*tman - 2*Bv0L*m3**2*m4**2*tman - Bv1L*m3**3*m4**2*tman - Bv2L*m3**3*m4**2*tman + 2*Bv3L*m3**3*m4**2*tman + Bs3L*m4**4*tman - 
     &        2*Bv0L*m4**4*tman - Bv1L*m3*m4**4*tman - Bv2L*m3*m4**4*tman + Bv3L*m3*m4**4*tman - Bs3L*m1**2*sman*tman - 2*Bv0L*m1**2*sman*tman - 
     &        Bs3L*m2**2*sman*tman - 2*Bv0L*m2**2*sman*tman + 2*Bv2R*m1*m2**2*sman*tman - 2*Bv3R*m1*m2**2*sman*tman + 2*Bs2R*m1*m3*sman*tman - 
     &        2*Bs3R*m1*m3*sman*tman + 4*Bv0R*m1*m3*sman*tman - Bv2L*m1**2*m3*sman*tman - Bv3L*m1**2*m3*sman*tman - Bv2L*m2**2*m3*sman*tman - 
     &        Bv3L*m2**2*m3*sman*tman - Bs3L*m3**2*sman*tman - 2*Bv0L*m3**2*sman*tman - Bv2L*m3**3*sman*tman - Bv3L*m3**3*sman*tman - Bs3L*m4**2*sman*tman - 
     &        2*Bv0L*m4**2*sman*tman - Bv2L*m3*m4**2*sman*tman - Bv3L*m3*m4**2*sman*tman + 2*Bv0L*sman**2*tman + Bv2L*m3*sman**2*tman + 2*Bv0L*m1**2*tman**2 + 
     &        4*Bv0L*m2**2*tman**2 + 2*Bv2R*m1*m2**2*tman**2 + 2*Bs2R*m1*m3*tman**2 + 4*Bv0R*m1*m3*tman**2 + Bv1L*m1**2*m3*tman**2 + Bv2L*m1**2*m3*tman**2 + 
     &        2*Bv2L*m2**2*m3*tman**2 - Bs3L*m3**2*tman**2 + 2*Bv0L*m3**2*tman**2 + Bv2L*m3**3*tman**2 - Bv3L*m3**3*tman**2 - Bs3L*m4**2*tman**2 + 
     &        4*Bv0L*m4**2*tman**2 + Bv1L*m3*m4**2*tman**2 + 2*Bv2L*m3*m4**2*tman**2 - Bv3L*m3*m4**2*tman**2 + Bs3L*sman*tman**2 + Bv3L*m3*sman*tman**2 - 
     &        2*Bv0L*tman**3 - Bv2L*m3*tman**3 + (-2*Bv2R*m1*m2**4 + Bv1L*m1**4*m3 - 4*Bs0L*m2**2*m3 - 2*Bs1R*m1*m2**2*m3 - 3*Bv1L*m1**2*m2**2*m3 + 
     &           Bv2L*m1**2*m2**2*m3 - Bv2L*m2**4*m3 + Bs3L*m1**2*m3**2 + Bs3L*m2**2*m3**2 + 2*Bv3R*m1*m2**2*m3**2 - 4*Bs0L*m3**3 - 2*Bs1R*m1*m3**3 + 
     &           3*Bv1L*m1**2*m3**3 - Bv3L*m1**2*m3**3 + Bv2L*m2**2*m3**3 + Bv3L*m2**2*m3**3 - Bs3L*m3**4 - Bv3L*m3**5 - 4*Bs0R*m1*m4**2 + Bs3L*m1**2*m4**2 + 
     &           Bs3L*m2**2*m4**2 - 4*Bv2R*m1*m2**2*m4**2 + 2*Bv3R*m1*m2**2*m4**2 - 8*Bs0L*m3*m4**2 - 4*Bs1R*m1*m3*m4**2 - 4*Bv0R*m1*m3*m4**2 + 
     &           2*Bv1L*m1**2*m3*m4**2 + Bv2L*m1**2*m3*m4**2 - Bv3L*m1**2*m3*m4**2 - Bv1L*m2**2*m3*m4**2 + Bv3L*m2**2*m3*m4**2 + 2*Bv3R*m1*m3**2*m4**2 + 
     &           Bv1L*m3**3*m4**2 + Bv2L*m3**3*m4**2 - 2*Bv3L*m3**3*m4**2 + Bs3L*m4**4 - 2*Bv2R*m1*m4**4 + 2*Bv3R*m1*m4**4 + Bv1L*m3*m4**4 + Bv2L*m3*m4**4 - 
     &           Bv3L*m3*m4**4 - 2*Bv1R*m1*(m2**4 + m4**4 + m1**2*(m2**2 + m4**2) + m2**2*(m3**2 + 3*m4**2 - sman)) - Bs3L*m1**2*sman - Bs3L*m2**2*sman - 
     &           2*Bv3R*m1*m2**2*sman + 4*Bs0L*m3*sman + 2*Bs1R*m1*m3*sman - Bv1L*m1**2*m3*sman + Bv3L*m1**2*m3*sman - Bv1L*m2**2*m3*sman - 
     &           Bv3L*m2**2*m3*sman + Bs3L*m3**2*sman - Bv1L*m3**3*sman + Bv3L*m3**3*sman - Bs3L*m4**2*sman - 2*Bv3R*m1*m4**2*sman - Bv1L*m3*m4**2*sman + 
     &           Bv3L*m3*m4**2*sman + Bv1L*m3*sman**2 + 2*Bv1R*m1*(m1**2 + m2**2 + m4**2)*tman + 
     &           (4*Bs0L*m3 + 2*Bs1R*m1*m3 + Bv1L*m1**2*m3 - Bv2L*m1**2*m3 + Bv1L*m2**2*m3 + Bv2L*m2**2*m3 - 2*Bv3R*m1*m3**2 + Bv1L*m3**3 - Bv2L*m3**3 - 
     &              2*Bv3R*m1*m4**2 + Bv1L*m3*m4**2 - Bv2L*m3*m4**2 + 4*Bv2R*m1*(m2**2 + m4**2) + 2*Bv3R*m1*sman)*tman - (2*Bv2R*m1 + Bv1L*m3)*tman**2 + 
     &           2*Bv0L*(m2**4 + m4**2*(m3**2 + m4**2 - 2*sman) + m1**2*(m2**2 + 3*m4**2 - tman) + m2**2*(3*m3**2 + 6*m4**2 - 2*sman - tman) - 
     &              (m3**2 + m4**2)*tman))*uman - (2*Bv1L*m1**2*m3 - Bv1L*m2**2*m3 + Bv2L*m2**2*m3 + Bs3L*m3**2 + Bv1L*m3**3 - Bv3L*m3**3 + Bs3L*m4**2 + 
     &           2*Bv1L*m3*m4**2 + Bv2L*m3*m4**2 - Bv3L*m3*m4**2 - Bs3L*sman + Bv3L*m3*sman + 2*Bv0L*(m2**2 + m4**2 - tman) - 
     &           2*Bv1R*m1*(m2**2 + m4**2 - tman) - Bv2L*m3*tman)*uman**2 + Bv1L*m3*uman**3 + 
     &        Bs2L*(m2**2 + m4**2 - tman)*((m2**2 + m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 - 2*m3**2 + sman - tman - uman) - 
     &           (m2**2 - m3**2 + m4**2)*uman + uman**2) + Bs1L*(m1**2 + m4**2 - uman)*
     &         ((m2**2 + m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 - 2*m3**2 + sman - tman - uman) - (m2**2 - m3**2 + m4**2)*uman + uman**2))
     &       *xv) + A3Rc*DRc*(-4*Bs1R*m1**2*m2**2 - 4*Bs2R*m1**2*m2**2 - 4*Bs3R*m1**2*m2**2 + 8*Bv0R*m1**2*m2**2 + 8*Bv1L*m1**3*m2**2 + 8*Bv0R*m2**4 + 
     &     4*Bv1L*m1*m2**4 - 4*Bv2L*m1*m2**4 - 4*Bv3L*m1*m2**4 - 16*Bs0R*m1**2*m3 + 8*Bs1L*m1**3*m3 - 8*Bs0R*m2**2*m3 + 4*Bs1L*m1*m2**2*m3 - 
     &     4*Bs2L*m1*m2**2*m3 - 4*Bs3L*m1*m2**2*m3 + 8*Bv0L*m1*m2**2*m3 - 4*Bv1R*m1**2*m2**2*m3 - 4*Bv2R*m1**2*m2**2*m3 - 4*Bv3R*m1**2*m2**2*m3 + 
     &     4*Bs1R*m1**2*m3**2 - 4*Bs2R*m1**2*m3**2 - 4*Bs3R*m1**2*m3**2 - 8*Bv0R*m1**2*m3**2 - 8*Bs2R*m2**2*m3**2 - 8*Bs3R*m2**2*m3**2 - 8*Bv0R*m2**2*m3**2 + 
     &     4*Bv1L*m1*m2**2*m3**2 + 4*Bv2L*m1*m2**2*m3**2 + 4*Bv3L*m1*m2**2*m3**2 - 8*Bs0R*m3**3 + 4*Bs1L*m1*m3**3 + 4*Bs2L*m1*m3**3 + 4*Bs3L*m1*m3**3 - 
     &     8*Bv0L*m1*m3**3 + 4*Bv1R*m1**2*m3**3 - 4*Bv2R*m1**2*m3**3 - 4*Bv3R*m1**2*m3**3 - 8*Bv2R*m2**2*m3**3 - 8*Bv3R*m2**2*m3**3 - 8*Bv0R*m2**2*sman - 
     &     4*Bv1L*m1*m2**2*sman + 8*Bs0R*m3*sman - 4*Bs1L*m1*m3*sman + 4*Bs2R*m3**2*sman + 4*Bs3R*m3**2*sman + 8*Bv0R*m3**2*sman + 4*Bv2R*m3**3*sman + 
     &     4*Bv3R*m3**3*sman + 4*Bs2R*m2**2*tman + 4*Bs3R*m2**2*tman - 4*Bv1L*m1*m2**2*tman + 8*Bs0R*m3*tman - 4*Bs1L*m1*m3*tman + 4*Bv2R*m2**2*m3*tman + 
     &     4*Bv3R*m2**2*m3*tman + 4*Bs1R*m1**2*uman + 4*Bs3R*m1**2*uman - 8*Bv0R*m1**2*uman - 4*Bv1L*m1**3*uman - 8*Bv0R*m2**2*uman - 4*Bv1L*m1*m2**2*uman + 
     &     4*Bv2L*m1*m2**2*uman + 8*Bv3L*m1*m2**2*uman - 4*Bs2L*m1*m3*uman + 4*Bs3L*m1*m3*uman + 8*Bv0L*m1*m3*uman - 4*Bv1R*m1**2*m3*uman + 
     &     4*Bv2R*m1**2*m3*uman + 4*Bv2R*m2**2*m3*uman + 4*Bs3R*m3**2*uman + 4*Bv3L*m1*m3**2*uman + 8*Bv0R*sman*uman + 4*Bv1L*m1*sman*uman - 
     &     4*Bv2R*m3*sman*uman - 4*Bs3R*tman*uman - 4*Bv3L*m1*uman**2 + 16*Bs0L*m1*(-m2**2 + m3**2 + uman) + 
     &     (2*Bv1L*m1**3*m2**4 + 2*Bv2L*m1*m2**6 + 4*Bs0R*m1**2*m2**2*m3 + 2*Bs1L*m1**3*m2**2*m3 + 2*Bv1R*m1**4*m2**2*m3 + 2*Bs2L*m1*m2**4*m3 + 
     &        4*Bv0L*m1*m2**4*m3 + 2*Bv2R*m1**2*m2**4*m3 - 2*Bs3R*m1**2*m2**2*m3**2 - 8*Bv0R*m1**2*m2**2*m3**2 + 2*Bv1L*m1**3*m2**2*m3**2 + 
     &        2*Bv2L*m1*m2**4*m3**2 - 2*Bv3L*m1*m2**4*m3**2 + 4*Bs0R*m1**2*m3**3 + 2*Bs1L*m1**3*m3**3 - 2*Bv1R*m1**4*m3**3 + 2*Bs2L*m1*m2**2*m3**3 - 
     &        2*Bs3L*m1*m2**2*m3**3 + 4*Bv0L*m1*m2**2*m3**3 - 2*Bv2R*m1**2*m2**2*m3**3 - 2*Bv3R*m1**2*m2**2*m3**3 + 2*Bs3R*m1**2*m3**4 - 
     &        2*Bv3L*m1*m2**2*m3**4 - 2*Bs3L*m1*m3**5 + 2*Bv3R*m1**2*m3**5 + 4*Bs0L*m1*m2**2*m4**2 - 2*Bs3R*m1**2*m2**2*m4**2 - 12*Bv0R*m1**2*m2**2*m4**2 + 
     &        4*Bv1L*m1**3*m2**2*m4**2 - 8*Bv0R*m2**4*m4**2 + 2*Bv1L*m1*m2**4*m4**2 + 6*Bv2L*m1*m2**4*m4**2 - 2*Bv3L*m1*m2**4*m4**2 + 16*Bs0R*m1**2*m3*m4**2 + 
     &        4*Bs1L*m1**3*m3*m4**2 + 8*Bs0R*m2**2*m3*m4**2 + 2*Bs1L*m1*m2**2*m3*m4**2 + 6*Bs2L*m1*m2**2*m3*m4**2 - 2*Bs3L*m1*m2**2*m3*m4**2 + 
     &        8*Bv0L*m1*m2**2*m3*m4**2 + 2*Bv1R*m1**2*m2**2*m3*m4**2 + 2*Bv2R*m1**2*m2**2*m3*m4**2 - 2*Bv3R*m1**2*m2**2*m3*m4**2 - 4*Bs0L*m1*m3**2*m4**2 + 
     &        2*Bs3R*m1**2*m3**2*m4**2 - 4*Bv0R*m1**2*m3**2*m4**2 - 8*Bv0R*m2**2*m3**2*m4**2 + 2*Bv1L*m1*m2**2*m3**2*m4**2 + 2*Bv2L*m1*m2**2*m3**2*m4**2 - 
     &        6*Bv3L*m1*m2**2*m3**2*m4**2 + 8*Bs0R*m3**3*m4**2 + 2*Bs1L*m1*m3**3*m4**2 + 2*Bs2L*m1*m3**3*m4**2 - 6*Bs3L*m1*m3**3*m4**2 + 
     &        8*Bv0L*m1*m3**3*m4**2 - 2*Bv1R*m1**2*m3**3*m4**2 - 2*Bv2R*m1**2*m3**3*m4**2 + 2*Bv3R*m1**2*m3**3*m4**2 - 8*Bv0R*m2**2*m4**4 + 
     &        4*Bv1L*m1*m2**2*m4**4 + 4*Bv2L*m1*m2**2*m4**4 - 4*Bv3L*m1*m2**2*m4**4 + 8*Bs0R*m3*m4**4 + 4*Bs1L*m1*m3*m4**4 + 4*Bs2L*m1*m3*m4**4 - 
     &        4*Bs3L*m1*m3*m4**4 + 8*Bv0L*m1*m3*m4**4 + 2*Bs3R*m1**2*m2**2*sman + 6*Bv0R*m1**2*m2**2*sman - 2*Bv1L*m1**3*m2**2*sman + 2*Bv0R*m2**4*sman - 
     &        2*Bv2L*m1*m2**4*sman + 2*Bv3L*m1*m2**4*sman - 4*Bs0R*m1**2*m3*sman - 2*Bs1L*m1**3*m3*sman + Bv1R*m1**4*m3*sman - 2*Bs2L*m1*m2**2*m3*sman + 
     &        2*Bs3L*m1*m2**2*m3*sman - 4*Bv0L*m1*m2**2*m3*sman + Bv1R*m1**2*m2**2*m3*sman + Bv2R*m1**2*m2**2*m3*sman + 2*Bv3R*m1**2*m2**2*m3*sman + 
     &        Bv2R*m2**4*m3*sman - 3*Bs3R*m1**2*m3**2*sman - Bs3R*m2**2*m3**2*sman + 2*Bv0R*m2**2*m3**2*sman + 4*Bv3L*m1*m2**2*m3**2*sman + 
     &        4*Bs3L*m1*m3**3*sman + Bv1R*m1**2*m3**3*sman - 3*Bv3R*m1**2*m3**3*sman + Bv2R*m2**2*m3**3*sman - Bv3R*m2**2*m3**3*sman - Bs3R*m3**4*sman - 
     &        Bv3R*m3**5*sman - Bs3R*m1**2*m4**2*sman + 2*Bv0R*m1**2*m4**2*sman - Bs3R*m2**2*m4**2*sman + 12*Bv0R*m2**2*m4**2*sman - 
     &        2*Bv1L*m1*m2**2*m4**2*sman - 2*Bv2L*m1*m2**2*m4**2*sman + 6*Bv3L*m1*m2**2*m4**2*sman - 8*Bs0R*m3*m4**2*sman - 2*Bs1L*m1*m3*m4**2*sman - 
     &        2*Bs2L*m1*m3*m4**2*sman + 6*Bs3L*m1*m3*m4**2*sman - 4*Bv0L*m1*m3*m4**2*sman + 2*Bv1R*m1**2*m3*m4**2*sman + Bv2R*m1**2*m3*m4**2*sman - 
     &        Bv3R*m1**2*m3*m4**2*sman + Bv1R*m2**2*m3*m4**2*sman + 2*Bv2R*m2**2*m3*m4**2*sman - Bv3R*m2**2*m3*m4**2*sman - 2*Bs3R*m3**2*m4**2*sman + 
     &        2*Bv0R*m3**2*m4**2*sman + Bv1R*m3**3*m4**2*sman + Bv2R*m3**3*m4**2*sman - 2*Bv3R*m3**3*m4**2*sman - Bs3R*m4**4*sman + 2*Bv0R*m4**4*sman + 
     &        Bv1R*m3*m4**4*sman + Bv2R*m3*m4**4*sman - Bv3R*m3*m4**4*sman + Bs3R*m1**2*sman**2 + Bs3R*m2**2*sman**2 - 2*Bv0R*m2**2*sman**2 - 
     &        2*Bv3L*m1*m2**2*sman**2 - 2*Bs3L*m1*m3*sman**2 - Bv1R*m1**2*m3*sman**2 + Bv3R*m1**2*m3*sman**2 - Bv2R*m2**2*m3*sman**2 + Bv3R*m2**2*m3*sman**2 + 
     &        2*Bs3R*m3**2*sman**2 + 2*Bv3R*m3**3*sman**2 + 2*Bs3R*m4**2*sman**2 - 2*Bv0R*m4**2*sman**2 - Bv1R*m3*m4**2*sman**2 - Bv2R*m3*m4**2*sman**2 + 
     &        2*Bv3R*m3*m4**2*sman**2 - Bs3R*sman**3 - Bv3R*m3*sman**3 - 2*Bv0R*m1**2*m2**2*tman - 2*Bv1L*m1**3*m2**2*tman - 2*Bv0R*m2**4*tman - 
     &        4*Bv2L*m1*m2**4*tman - 4*Bs0R*m1**2*m3*tman - 2*Bs1L*m1**3*m3*tman - Bv1R*m1**4*m3*tman - 4*Bs2L*m1*m2**2*m3*tman - 8*Bv0L*m1*m2**2*m3*tman - 
     &        Bv1R*m1**2*m2**2*m3*tman - 3*Bv2R*m1**2*m2**2*m3*tman - Bv2R*m2**4*m3*tman + Bs3R*m1**2*m3**2*tman + 4*Bv0R*m1**2*m3**2*tman + 
     &        Bs3R*m2**2*m3**2*tman - 2*Bv0R*m2**2*m3**2*tman - 2*Bv2L*m1*m2**2*m3**2*tman + 2*Bv3L*m1*m2**2*m3**2*tman - 2*Bs2L*m1*m3**3*tman + 
     &        2*Bs3L*m1*m3**3*tman - 4*Bv0L*m1*m3**3*tman - Bv1R*m1**2*m3**3*tman + 2*Bv2R*m1**2*m3**3*tman + Bv3R*m1**2*m3**3*tman - Bv2R*m2**2*m3**3*tman + 
     &        Bv3R*m2**2*m3**3*tman + Bs3R*m3**4*tman + Bv3R*m3**5*tman + Bs3R*m1**2*m4**2*tman - 2*Bv0R*m1**2*m4**2*tman + Bs3R*m2**2*m4**2*tman + 
     &        4*Bv0R*m2**2*m4**2*tman - 2*Bv1L*m1*m2**2*m4**2*tman - 6*Bv2L*m1*m2**2*m4**2*tman + 2*Bv3L*m1*m2**2*m4**2*tman - 8*Bs0R*m3*m4**2*tman - 
     &        2*Bs1L*m1*m3*m4**2*tman - 6*Bs2L*m1*m3*m4**2*tman + 2*Bs3L*m1*m3*m4**2*tman - 12*Bv0L*m1*m3*m4**2*tman - 2*Bv1R*m1**2*m3*m4**2*tman - 
     &        Bv2R*m1**2*m3*m4**2*tman + Bv3R*m1**2*m3*m4**2*tman - Bv1R*m2**2*m3*m4**2*tman - 2*Bv2R*m2**2*m3*m4**2*tman + Bv3R*m2**2*m3*m4**2*tman + 
     &        2*Bs3R*m3**2*m4**2*tman - 2*Bv0R*m3**2*m4**2*tman - Bv1R*m3**3*m4**2*tman - Bv2R*m3**3*m4**2*tman + 2*Bv3R*m3**3*m4**2*tman + Bs3R*m4**4*tman - 
     &        2*Bv0R*m4**4*tman - Bv1R*m3*m4**4*tman - Bv2R*m3*m4**4*tman + Bv3R*m3*m4**4*tman - Bs3R*m1**2*sman*tman - 2*Bv0R*m1**2*sman*tman - 
     &        Bs3R*m2**2*sman*tman - 2*Bv0R*m2**2*sman*tman + 2*Bv2L*m1*m2**2*sman*tman - 2*Bv3L*m1*m2**2*sman*tman + 2*Bs2L*m1*m3*sman*tman - 
     &        2*Bs3L*m1*m3*sman*tman + 4*Bv0L*m1*m3*sman*tman - Bv2R*m1**2*m3*sman*tman - Bv3R*m1**2*m3*sman*tman - Bv2R*m2**2*m3*sman*tman - 
     &        Bv3R*m2**2*m3*sman*tman - Bs3R*m3**2*sman*tman - 2*Bv0R*m3**2*sman*tman - Bv2R*m3**3*sman*tman - Bv3R*m3**3*sman*tman - Bs3R*m4**2*sman*tman - 
     &        2*Bv0R*m4**2*sman*tman - Bv2R*m3*m4**2*sman*tman - Bv3R*m3*m4**2*sman*tman + 2*Bv0R*sman**2*tman + Bv2R*m3*sman**2*tman + 2*Bv0R*m1**2*tman**2 + 
     &        4*Bv0R*m2**2*tman**2 + 2*Bv2L*m1*m2**2*tman**2 + 2*Bs2L*m1*m3*tman**2 + 4*Bv0L*m1*m3*tman**2 + Bv1R*m1**2*m3*tman**2 + Bv2R*m1**2*m3*tman**2 + 
     &        2*Bv2R*m2**2*m3*tman**2 - Bs3R*m3**2*tman**2 + 2*Bv0R*m3**2*tman**2 + Bv2R*m3**3*tman**2 - Bv3R*m3**3*tman**2 - Bs3R*m4**2*tman**2 + 
     &        4*Bv0R*m4**2*tman**2 + Bv1R*m3*m4**2*tman**2 + 2*Bv2R*m3*m4**2*tman**2 - Bv3R*m3*m4**2*tman**2 + Bs3R*sman*tman**2 + Bv3R*m3*sman*tman**2 - 
     &        2*Bv0R*tman**3 - Bv2R*m3*tman**3 + (-2*Bv2L*m1*m2**4 + Bv1R*m1**4*m3 - 4*Bs0R*m2**2*m3 - 2*Bs1L*m1*m2**2*m3 - 3*Bv1R*m1**2*m2**2*m3 + 
     &           Bv2R*m1**2*m2**2*m3 - Bv2R*m2**4*m3 + Bs3R*m1**2*m3**2 + Bs3R*m2**2*m3**2 + 2*Bv3L*m1*m2**2*m3**2 - 4*Bs0R*m3**3 - 2*Bs1L*m1*m3**3 + 
     &           3*Bv1R*m1**2*m3**3 - Bv3R*m1**2*m3**3 + Bv2R*m2**2*m3**3 + Bv3R*m2**2*m3**3 - Bs3R*m3**4 - Bv3R*m3**5 - 4*Bs0L*m1*m4**2 + Bs3R*m1**2*m4**2 + 
     &           Bs3R*m2**2*m4**2 - 4*Bv2L*m1*m2**2*m4**2 + 2*Bv3L*m1*m2**2*m4**2 - 8*Bs0R*m3*m4**2 - 4*Bs1L*m1*m3*m4**2 - 4*Bv0L*m1*m3*m4**2 + 
     &           2*Bv1R*m1**2*m3*m4**2 + Bv2R*m1**2*m3*m4**2 - Bv3R*m1**2*m3*m4**2 - Bv1R*m2**2*m3*m4**2 + Bv3R*m2**2*m3*m4**2 + 2*Bv3L*m1*m3**2*m4**2 + 
     &           Bv1R*m3**3*m4**2 + Bv2R*m3**3*m4**2 - 2*Bv3R*m3**3*m4**2 + Bs3R*m4**4 - 2*Bv2L*m1*m4**4 + 2*Bv3L*m1*m4**4 + Bv1R*m3*m4**4 + Bv2R*m3*m4**4 - 
     &           Bv3R*m3*m4**4 - 2*Bv1L*m1*(m2**4 + m4**4 + m1**2*(m2**2 + m4**2) + m2**2*(m3**2 + 3*m4**2 - sman)) - Bs3R*m1**2*sman - Bs3R*m2**2*sman - 
     &           2*Bv3L*m1*m2**2*sman + 4*Bs0R*m3*sman + 2*Bs1L*m1*m3*sman - Bv1R*m1**2*m3*sman + Bv3R*m1**2*m3*sman - Bv1R*m2**2*m3*sman - 
     &           Bv3R*m2**2*m3*sman + Bs3R*m3**2*sman - Bv1R*m3**3*sman + Bv3R*m3**3*sman - Bs3R*m4**2*sman - 2*Bv3L*m1*m4**2*sman - Bv1R*m3*m4**2*sman + 
     &           Bv3R*m3*m4**2*sman + Bv1R*m3*sman**2 + 2*Bv1L*m1*(m1**2 + m2**2 + m4**2)*tman + 
     &           (4*Bs0R*m3 + 2*Bs1L*m1*m3 + Bv1R*m1**2*m3 - Bv2R*m1**2*m3 + Bv1R*m2**2*m3 + Bv2R*m2**2*m3 - 2*Bv3L*m1*m3**2 + Bv1R*m3**3 - Bv2R*m3**3 - 
     &              2*Bv3L*m1*m4**2 + Bv1R*m3*m4**2 - Bv2R*m3*m4**2 + 4*Bv2L*m1*(m2**2 + m4**2) + 2*Bv3L*m1*sman)*tman - (2*Bv2L*m1 + Bv1R*m3)*tman**2 + 
     &           2*Bv0R*(m2**4 + m4**2*(m3**2 + m4**2 - 2*sman) + m1**2*(m2**2 + 3*m4**2 - tman) + m2**2*(3*m3**2 + 6*m4**2 - 2*sman - tman) - 
     &              (m3**2 + m4**2)*tman))*uman - (2*Bv1R*m1**2*m3 - Bv1R*m2**2*m3 + Bv2R*m2**2*m3 + Bs3R*m3**2 + Bv1R*m3**3 - Bv3R*m3**3 + Bs3R*m4**2 + 
     &           2*Bv1R*m3*m4**2 + Bv2R*m3*m4**2 - Bv3R*m3*m4**2 - Bs3R*sman + Bv3R*m3*sman + 2*Bv0R*(m2**2 + m4**2 - tman) - 
     &           2*Bv1L*m1*(m2**2 + m4**2 - tman) - Bv2R*m3*tman)*uman**2 + Bv1R*m3*uman**3 + 
     &        Bs2R*(m2**2 + m4**2 - tman)*((m2**2 + m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 - 2*m3**2 + sman - tman - uman) - 
     &           (m2**2 - m3**2 + m4**2)*uman + uman**2) + Bs1R*(m1**2 + m4**2 - uman)*
     &         ((m2**2 + m3**2 + m4**2 - sman - tman)*(sman - tman) + m1**2*(2*m2**2 - 2*m3**2 + sman - tman - uman) - (m2**2 - m3**2 + m4**2)*uman + uman**2))
     &       *xv) ))/(uman - mq3c**2))
                  

  
	end

