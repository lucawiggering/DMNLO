#define DMNLO_VERSION "1.0"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <map>
#include <getopt.h>
#include <stdbool.h>


using namespace std;



struct Args {
    char *logFile;
    string inputFile, slhaFile;
    double pcm,muR;
    long na,nb,n1,n2; // process a + b --> 1 + 2
    int fluxFlag,renscheme,stopAfterLO = 0,stopAfterNLO = 0,SommerfeldAlone=0, full = 0,legacy=0, choosesol = 0,formfactor = 0;
    bool correct_ChiChi2QQ = true,correct_NeuQ2qx = true, correct_staustau2QQ = true,
         correct_stst2QQ = true, correct_QQ2xx = true, correct_stsT2xx = true;
    bool useDD = false;
    map<string, string> arguments,config;
};

extern void printBanner();
extern map<string, string> readInputFile(string filename);
extern void parseArgs(int argc, char *argv[], Args &args);
extern double compute(Args &args, int verbosity);
void computeDD(Args &args, double *dd_results);


// DM@NLO processes
extern "C"
{
    void chichi2qq_(double &renscale, double &pcm, int *mflags, double &cut, double *result, const char *slhafilename);
    void neuq2qx_(double &renscale, long &n1,long &n2,long &n3,long &n4, double &pcm, int *mflags, double *result, double *error, const char *slhafilename);
    void staustau2qq_(double &renscale,double &pcm,int *mflags,double *result,const char *slhafilename);
    void stst2qq_(double &renscale,double &pcm,int *mflags,double *result, double *error, const char *slhafilename);
    void qq2xx_(double &renscale,long &n1,long &n2,long &n3,long &n4,long &n5,double &pcm,int *mflags,double *result, double *error,const char *slhafilename);
    void stst2xx_(double &renscale, double &pcm, int *mflags, double *result, double *error, double &gaugek1, double &gaugek2,const char *slhafilename);
    void directdetection_(double &pcm, double &costh, int &inLoop_Mode, int &inDenom_Mode, int &inBox_Mode, int &inNucl_Mode, double &Qout, double *results, int *iflags, const char *slhafilename);

    int check_renscheme_(int &ischeme, double &renscale, const char *slhafilename);
}