/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about SUSY and Higgs masses 
      */      


#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */      
        
#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  


/*===== end of Modules  ======*/


/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


//  Example file for using DM@NLO. File should be disposed in MSSM/lib

// --- input parameters for DM@NLO ---
static double muR; // renormalization scale
static char *slhaFile; // SLHA file
static int renscheme; // renormalization scheme
static int RDoption; // This defines the relic density (RD) calculation option ( 0 = MO result (CalcHEP)
//                                                                               1 = LO (DM@NLO) 
//                                                                               2 = NLO (DM@NLO) 
//                                                                               3 = full (NLO + Sommerfeld if Sommerfeld available) (DM@NLO) 
//                                                                               4 = Sommerfeld alone (DM@NLO) )
static double ChannelFilter = 0.02; // only (co)annihlation channels above the value given by ChannelFiler are corrected by DM@NLO



extern double cs_dmnlo(int order,long na, long nb, long n1, long n2, double pcm, double muR, char *slha, int renscheme, int choosesol, int *corrFlags);
extern void dd_dmnlo(int order, double muR, char *slha, int renscheme, int choosesol, int formfactor, double *cs); // for direct detection, please cite 1607.06396
extern int canImprove_dmnlo(long n1, long n2, long n3, long n4);
extern int consistent_RS_dmnlo(int renscheme, char *slha, double muR);

// function that computes the tree-level stop masses in order to compute QSUSY
// which we define as the geomeetric mean of the two stop masses
double computeQSUSY(){
   double mZ,sW,mt,atop,mu,tb,mqL3,mtR,QSUSY;
   mZ = findValW("MZ");
   sW = findValW("SW");
   mt = findValW("Mtp");

   atop = findValW("At");
   mu = findValW("mu");
   tb = findValW("tb");
   mqL3 = findValW("Mq3");
   mtR = findValW("Mu3");

   double MassStop[2][2],trace,det,mStopTree1,mStopTree2;
   MassStop[0][0] = pow(mqL3,2)+pow(mZ,2)*(0.5-2./3.*pow(sW,2))*cos(2*atan(tb))+pow(mt,2); //M_LL
   MassStop[1][1] = pow(mtR,2)+2./3.*pow(mZ,2)*cos(2*atan(tb))*pow(sW,2)+pow(mt,2); //M_RR
   MassStop[0][1] = -mt*(atop+mu/tb); //M_LR
   MassStop[1][0] = -mt*(atop+mu/tb); //M_RL

   trace = MassStop[0][0]+MassStop[1][1];
   det = MassStop[0][0]*MassStop[1][1]-MassStop[0][1]*MassStop[1][0];
   mStopTree1 = sqrt(0.5*(trace -sqrt(pow(trace,2)-4*det)));
   mStopTree2 = sqrt(0.5*(trace +sqrt(pow(trace,2)-4*det)));

   QSUSY = sqrt(mStopTree1*mStopTree2);
   return QSUSY; 
}

// function that checks whether the contribution from the (co)annihilation process n1+n2 --> n3+n4 is larger than the value specified in ChannelFilter
// this boundary is set in this example to 2%
bool shallImprove(long n1, long n2, long n3, long n4)
{
   for(int i=0; omegaCh[i].weight>ChannelFilter  ;i++)
   {												
      if (strcmp(omegaCh[i].prtcl[0], pdg2name(n1))==0 && strcmp(omegaCh[i].prtcl[1], pdg2name(n2))==0
         && strcmp(omegaCh[i].prtcl[2], pdg2name(n3))==0 && strcmp(omegaCh[i].prtcl[3], pdg2name(n4))==0) {
         return true;
      }
   }
   return false;
}

void improveCrossSection(long n1,long n2,long n3,long n4 ,double PcmIn,  double * res) 
{
   if(!RDoption) return;
   if(!shallImprove(n1, n2, n3, n4)) return;
   
   int order,corrFlag[6],choosesol = 0; 
   double result=0.;
   // --- Specify the order of the calculation ---
   order = RDoption-1;

   // corrFlag allows to decide what kind of processes are replaced by DM@NLO
   // 0 = don't replace these processes, 1 = replace these processes
   corrFlag[0] = 1; // ChiChi2QQ:  gaugino (co-)annihilation into quarks, please cite 0709.0043, 0901.0481, 0907.0030, 1404.2931
   corrFlag[1] = 1; // NeuQ2qx: neutralino-stop coannihilation, please cite 1212.5241, 1409.2898
   corrFlag[2] = 1; // staustau2QQ: stau annihilation into top quarks, please cite 1909.09527
   corrFlag[3] = 1; // stst2QQ: squark-pair annihilation into heavy quarks, please cite 1903.10998 
   corrFlag[4] = 1; // QQ2xx: stop annihilation into electroweak final states, please cite 1410.8063
   corrFlag[5] = 1; // stsT2xx: stop-antistop annihilation into gluons and light quarks, please cite 2210.05260

   // Calculate cross section 
   result = cs_dmnlo(order,n1,n2,n3,n4,PcmIn,muR,slhaFile,renscheme,choosesol,corrFlag);

   if(result > 0) *res = result;  

   return;
}



int main(int argc,char** argv){  
   int err;
   char cdmName[10];
   int spin2, charge3,cdim;

// sysTimeLim=1000; 
   ForceUG=0;   /* to Force Unitary Gauge assign 1 */
   //useSLHAwidth=0;
//  nPROCSS=0; /* to switch off multiprocessor calculations */

   printf("\n========= SLHA file input =========\n");

   if(argc <2) 
   {  
      printf("The program needs one argument:the name of SLHA input file.\n"
            "Example: ./main suspect2_lha.out \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
   err=lesHinput(argv[1]);
   if(err) exit(2);

         
   if(err==-1)     { printf("Can not open the file\n"); exit(2);}
   else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}

   err=sortOddParticles(cdmName);

   if(err) { printf("Can't calculate %s\n",cdmName); return 1;}


   VWdecay=1; VZdecay=1;


   qNumbers(cdmName,&spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2  mass=%.2E\n",
   cdmName,       spin2, Mcdm); 

   if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
   if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                              else o1Contents(stdout);

                             
#ifdef MASSES_INFO
   printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
   printHiggs(stdout);
   printMasses(stdout,1);  
#endif

// ----- DM@NLO settings -------------
   muR = computeQSUSY(); // renormalization scale 
   renscheme = 1; // renormalization scheme, for more info see e.g. the settings file DMNLO.in or the manual
   slhaFile = argv[1]; // SLHA file

   int choosesol = 0, formfactor = 0;

   int renscheme_works = consistent_RS_dmnlo(renscheme,slhaFile,muR);
// ------------------------------------

#ifdef OMEGA
   int fast=0,i;
   double Beps=1.E-5, cut=0.01,w,w_dmnlo=0.;
   double Omega,Xf=25;

// --- CalcHEP result ---
   printf("\n==== Calculation of relic density with CalcHEP =====\n"); 
   RDoption = 0;
   Omega=darkOmega(&Xf,fast,Beps,&err);                                       
   printf("CalcHEP: Xf=%.2e Omega=%.3e err=%d\n",Xf,Omega,err);
   if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);

   // print percentage of (co)annihilation channels that can be corrected with dmnlo
    w = 1.;
    i = 0;
    while (w>ChannelFilter) {
        if(canImprove_dmnlo(pNum(omegaCh[i].prtcl[0]),pNum(omegaCh[i].prtcl[1]),pNum(omegaCh[i].prtcl[2]),pNum(omegaCh[i].prtcl[3])))  w_dmnlo += omegaCh[i].weight;
        i++;
        w = omegaCh[i].weight;
    }
    printf("%d%% of the contributions can be corrected with DM@NLO\n",(int) round(w_dmnlo*100));

// --- DM@NLO tree-level result ---
   printf("\n==== Calculation of relic density with DM@NLO @ LO =====\n"); 
   RDoption = 1; Omega = -999;
   if(renscheme_works) Omega=darkOmega(&Xf,fast,Beps,&err);                                       
   printf("DM@NLO tree: Xf=%.2e Omega=%.3e err=%d\n",Xf,Omega,err);
   if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);

// --- DM@NLO full result ---
   printf("\n==== Calculation of relic density with DM@NLO @ NLO and beyond =====\n"); 
   RDoption = 3; Omega = -999;
   if(renscheme_works) Omega=darkOmega(&Xf,fast,Beps,&err);                                       
   printf("DM@NLO tree: Xf=%.2e Omega=%.3e err=%d\n",Xf,Omega,err);
   if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);

   VZdecay=1; VWdecay=1; cleanDecayTable();
#endif
 


#ifdef CDM_NUCLEON
   int order;
   double ddres[4],pval;
   char *expName="";
/// --- MicrOMEGAs ---
   printf("\n==== Calculation of CDM-nucleons amplitudes with MicrOMEGAs  =====\n");   
   double pA0[2],pA5[2],nA0[2],nA5[2];
   double Nmass=0.939; /*nucleon mass*/
   double SCcoeff;        
   double csSIp,csSIn,csSDp,csSDn;
   double pb2cm2 = 1e-36; // convert pb to cm^2
   int sI,sD; 

   nucleonAmplitudes(CDM1,pA0,pA5,nA0,nA5);
   SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
   csSIp=  SCcoeff*pA0[0]*pA0[0];  
   csSDp=3*SCcoeff*pA5[0]*pA5[0];  
   csSIn=  SCcoeff*nA0[0]*nA0[0];  
   csSDn=3*SCcoeff*nA5[0]*nA5[0];  
                    
   printf("\n==== %s-nucleon cross sections [cm^2] ====\n",CDM1);
   printf(" proton  SI %.3E  SD %.3E\n",csSIp*pb2cm2,csSDp*pb2cm2);
   printf(" neutron SI %.3E  SD %.3E\n",csSIn*pb2cm2,csSDn*pb2cm2);

   if(pA0[0]*nA0[0]<0) sI=-1; else sI=1;
   if(pA5[0]*nA5[0]<0) sD=-1; else sD=1;

   pval=DD_pvalCS(AllDDexp, Maxwell, csSIp, sI*csSIn,csSDp,sD*csSDn, &expName);
   if(pval<0.1 )  printf("MicrOMEGAs: Excluded by %s [CDM_NUCLEON] %.1f%% \n", expName, 100*(1-pval)); 
   else printf("MicrOMEGAs: Not excluded by DD experiments  at 90%% level \n");        

/// --- DM@NLO tree-level ---
   order = 0; for(int i=0;i<4;i++) ddres[i] = -999;
   printf("\n==== Calculation of CDM-nucleons amplitudes with DM@NLO @ LO  =====\n");  
   if(renscheme_works) dd_dmnlo(order, muR, slhaFile, renscheme, choosesol, formfactor, ddres);   
   printf("\n==== %s-nucleon cross sections [cm^2] ====\n",CDM1);
   printf(" proton  SI %.3E  SD %.3E\n",ddres[0],ddres[1]);
   printf(" neutron SI %.3E  SD %.3E\n",ddres[2],ddres[3]);

   pval=DD_pvalCS(AllDDexp, Maxwell, csSIp, sI*csSIn,csSDp,sD*csSDn, &expName);
   if(pval<0.1 )  printf("DM@NLO Tree: Excluded by %s [CDM_NUCLEON] %.1f%% \n", expName, 100*(1-pval)); 
   else printf("DM@NLO Tree: Not excluded by DD experiments  at 90%% level \n");        

/// --- DM@NLO NLO ---
   order = 1; for(int i=0;i<4;i++) ddres[i] = -999;
   printf("\n==== Calculation of CDM-nucleons amplitudes with DM@NLO @ NLO  =====\n"); 
   if(renscheme_works) dd_dmnlo(order, muR, slhaFile, renscheme, choosesol, formfactor, ddres);    
   printf("\n==== %s-nucleon cross sections [cm^2] ====\n",CDM1);
   printf(" proton  SI %.3E  SD %.3E\n",ddres[0],ddres[1]);
   printf(" neutron SI %.3E  SD %.3E\n",ddres[2],ddres[3]);

   pval=DD_pvalCS(AllDDexp, Maxwell, csSIp, sI*csSIn,csSDp,sD*csSDn, &expName);
   if(pval<0.1 )  printf("DM@NLO NLO: Excluded by %s [CDM_NUCLEON] %.1f%% \n", expName, 100*(1-pval)); 
   else printf("DM@NLO NLO: Not excluded by DD experiments  at 90%% level \n");  
        
#endif
  


   return 0;
}