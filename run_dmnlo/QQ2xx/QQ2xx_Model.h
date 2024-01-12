* model.h

#include "../util/DMNLO_Model.h"

* Additional Parameters	

        double precision zerocheck

         parameter (zerocheck = 1.D-30)


* Additional code Indizes and masses

	integer sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2,f3type,f4type,f3gen,f4gen,xtype,ytype,ztype,xhtype,yhtype, 
     &          propgent,propgenu,propgentc,propgenuc,proptypet,proptypeu,proptypetc,proptypeuc,Looptype0,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2,
     &          Loopmeigent,Looptypet,Loopgent,Loopmeigenu,Looptypeu,Loopgenu,sumstartH,sumendH,sumstartV,sumstartVV,sumendV,n1,n2,n3,n4,n5,
     &          channel,corrProp,corrVert,corrBox,isomm

	double precision Pcm,m1,m2,m3,m4,m5,M0Loop,M1Loop,M2Loop
	double precision ef_,I3f_,hf_
	double precision testvar,testvar1,testvar2,testvar3,testvar4
	double precision treeperc(7)
        double precision MV(4),MNeuCha2(6)
         
        common /add_para/ m1,m2,m3,m4,m5,M0Loop,M1Loop,M2Loop
	common /add_para/ testvar,testvar1,testvar2,testvar3,testvar4
        common /add_para/ treeperc
        common /add_para/ sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2,f3type,f4type,f3gen,f4gen,xtype,ytype,ztype,xhtype,yhtype,
     &                    propgent,propgenu,propgentc,propgenuc,proptypet,proptypeu,proptypetc,proptypeuc,Looptype0,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2,
     &                    Loopmeigent,Looptypet,Loopgent,Loopmeigenu,Looptypeu,Loopgenu,sumstartH,sumendH,sumstartV,sumstartVV,sumendV
     
        common /sm_para/ MV
        common /mssm_para/ MNeuCha2
        common /global/ Pcm,n1,n2,n3,n4,n5,channel,corrProp,corrVert,corrBox
        common /sommerfeld/ isomm
