       
        integer Looptypeout1,Looptypeout2,Loopgenout1,Loopgenout2
        double precision sfm,sfmc,colorfactorSSVV,colorfactorSSVGl,colorfactorSSSV,colorfactorSSSSV,colorfactorSSVVV,colorfactorSSSGl,colorfactorQQQQ,colorfactorQQVV,
     &                   colorfactorVS,colorfactorV,colorfactorS,colorfactorGl,colorfactorTU,colorfactorQ,colorfactorSH,colorfactorSV,colorfactorVH_HV,
     &                   colorfactorVV_SS,colorfactorVV_TT,colorfactorVV_TU,colorfactorVV_QQ,colorfactorVV_QTU,colorfactorVV_QS
        complex*16 DenSmH(6),DenSmV(4),Coupl_A,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_APc,Coupl_AMc,Coupl_B,Coupl_BP,Coupl_BM,Coupl_Bc,Coupl_BPc,Coupl_BMc,propa,propac,SS,colorfactorVV_STU        
        complex*16 SST(7,2),SSU(7,2),Coupl_AQS(5),Coupl_AQV(5),Coupl_AS(6,6),Coupl_APV(8,4),Coupl_AMV(8,4),Coupl_AT(6,2),Coupl_BT(6,2),
     &             Coupl_APT(8,2),Coupl_BPT(8,2),Coupl_AMT(8,2),Coupl_BMT(8,2),Coupl_AU(6,2),Coupl_BU(6,2),Coupl_APU(8,2),Coupl_BPU(8,2),Coupl_AMU(8,2),Coupl_BMU(8,2)
         
        common /Init_2to2/    Looptypeout1,Looptypeout2,Loopgenout1,Loopgenout2,DenSmH,DenSmV,Coupl_A,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_APc,
     &                        Coupl_AMc,Coupl_B,Coupl_BP,Coupl_BM,Coupl_Bc,Coupl_BPc,Coupl_BMc,propa,propac,sfm,sfmc,colorfactorSSVV,colorfactorSSVGl,colorfactorSSSV,colorfactorSSSSV,
     &                        colorfactorSSVVV,colorfactorSSSGl,colorfactorQQQQ,colorfactorQQVV,colorfactorVS,colorfactorV,colorfactorS,colorfactorGl,colorfactorTU,
     &                        colorfactorQ,colorfactorSH,colorfactorSV,colorfactorVH_HV,colorfactorVV_SS,colorfactorVV_TT,colorfactorVV_TU,colorfactorVV_QQ,colorfactorVV_QTU,
     &                        colorfactorVV_QS,colorfactorVV_STU,SS,SST,SSU,Coupl_AQS,Coupl_AQV,Coupl_AS,Coupl_APV,Coupl_AMV,Coupl_AT,Coupl_BT,
     &                        Coupl_APT,Coupl_BPT,Coupl_AMT,Coupl_BMT,Coupl_AU,Coupl_BU,Coupl_APU,Coupl_BPU,Coupl_AMU,Coupl_BMU 


