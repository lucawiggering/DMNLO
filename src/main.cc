#include "DMNLO.h"



int main(int argc, char *argv[]) {
  printBanner();
  struct Args args;
  parseArgs(argc, argv, args);
  if(args.useDD){
    double r[4];
    computeDD(args,r);

    printf("--- spin-independent dark matter-nucleon scattering cross section in cm^2 ---\n");
    printf("proton SI : %e   neutron SI : %e\n",r[0],r[2]);
    printf("--- spin-dependent dark matter-nucleon cross sectino in cm^2 ---\n");
    printf("proton SD : %e   neutron SD : %e\n",r[1],r[3]);
  }
  else{
    double r = compute(args,0);
    printf("Result[GeV^-2] = %e\n",r);
  }

}

