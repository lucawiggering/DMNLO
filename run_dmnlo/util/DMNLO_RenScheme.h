        double precision dZG, DeltaGs, dZC
        double complex dZfL(4,3), dZfR(4,3), dZSf(2,2,4,3), dMf(4,3), dMSf2(2,4,3), dUSf(2,2,4,3), dAf(4,3),dZGlL,dZGlR,dMGl
        double complex PifVL(4,3), PifVR(4,3), PifSL(4,3), PifSR(4,3), dPifVL(4,3), dPifVR(4,3), dPifSL(4,3), dPifSR(4,3)
        double complex PiglVL, PiglVR, PiglSL, PiglSR, dPiglVL, dPiglVR, dPiglSL, dPiglSR
        double complex PifVLth(4,3), PifVRth(4,3), PifSLth(4,3), PifSRth(4,3), dMfth(4,3)
        double complex PiSf(2,2,4,3), dPiSf(2,2,4,3)
        double complex PifVLtest(4,3), PifVRtest(4,3), PifSLtest(4,3), PifSRtest(4,3)
        double complex dPiGT,dPiC
        double complex choosesol
        integer renfail

* clean up the different bottom masses	        
        double precision DeltaQCD, Iabc, DeltaMb, DeltaMbQCD 
        integer Res_Style, NNLO_Res

        common /counterterm/ dZfL,dZfR,dZSf,dZG,dZC,dMf,dMSf2,dUSf,dAf,DeltaGs,dZGlL,dZGlR,dMGl
        common /fermionPi/ PifVL, PifVR, PifSL, PifSR, dPifVL, dPifVR, dPifSL, dPifSR
        common /fermionPi/ PiglVL, PiglVR, PiglSL, PiglSR, dPiglVL, dPiglVR, dPiglSL, dPiglSR
        common /fermionPith/ PifVLth, PifVRth, PifSLth, PifSRth, dMfth
        common /sfermionPi/ PiSf, dPiSf
        common /gluonPi/ dPiGT,dPiC
        common /HOcorr/ DeltaQCD, Iabc, DeltaMb, DeltaMbQCD

        common /fermionPitest/ PifVLtest, PifVRtest, PifSLtest, PifSRtest
        common /DelMb/ Res_Style, NNLO_Res
        common /USfchoose/ choosesol
        common /renfail/ renfail
