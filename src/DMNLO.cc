#include "DMNLO.h"
#include <iostream>

// Gets option from command-line if set, from input file otherwise.
#define get_option(name) (args.arguments[name] == "" ? args.config[name] : args.arguments[name])


void printBanner(){
   printf("\n"
         "  Welcome to DMNLO " DMNLO_VERSION "\n"
         "  Copyright 2007-2023 Björn Herrmann. \n"
         "  Copyright 2007-2023 Karol Kovarik. \n"
         "  Copyright 2020-2023 Luca Wiggering. \n"
         "\n"
         "  Licensed under the EUPL 1.1 or later.\n"
         "  See https://dmnlo.hepforge.org for more information.\n"
         "\n"
          "  For general DM@NLO usage, please cite:\n"
         "  - J. Harz, B. Herrmann, M. Klasen, K. Kovarik and L. P. Wiggering [2312.17206].\n"
         "\n"
         "  For direct detection, please cite:\n"
         "  - M. Klasen, K. Kovarik and P. Steppeler, Phys. Rev. D 94 (2016) 095002, [1607.06396]\n"
         "\n"
         "  For gaugino (co-)annihilation into quarks, please cite:\n"
         "  - B. Herrmann and M. Klasen, Phys. Rev. D 76 (2007) 117704, [0709.0043].\n"
         "  - B. Herrmann, M. Klasen and K. Kovarik, Phys. Rev. D 79 (2009) 061701, [0901.0481].\n"
         "  - B. Herrmann, M. Klasen and K. Kovarik, Phys. Rev. D 80 (2009) 085025, [0907.0030].\n"
         "  - B. Herrmann, M. Klasen, K. Kovarik, M. Meinecke and P. Steppeler,  Phys. Rev. D 89 (2014) 114012, [1404.2931].\n"
         "\n"
         "  For neutralino-stop coannihilation, please cite:\n"
         "  - J. Harz, B. Herrmann, M. Klasen, K. Kovarik and Q. L.Boulc’h, Phys. Rev. D 87 (2013) 054031, [1212.5241].\n"
         "  - J. Harz, B. Herrmann, M. Klasen and K. Kovarik, Phys. Rev. D 91 (2015) 034028, [1409.2898].\n"
         "\n"
         "  For stau annihilation into top quarks, please cite:\n"
         "  - J. Branahl, J. Harz, B. Herrmann, M. Klasen, K. Kovařík and S. Schmiemann, Phys. Rev. D 100 (2019) 115003, [1909.09527].\n"
         "\n"
         "  For squark-pair annihilation into heavy quarks, please cite:\n"
         "  - S. Schmiemann, J. Harz, B. Herrmann, M. Klasen and K. Kovařík, Phys. Rev. D 99 (2019) 095015, [1903.10998].\n"
         "\n"
         "  For stop annihilation into electroweak final states, please cite:\n"
         "  - J. Harz, B. Herrmann, M. Klasen, K. Kovařík and M. Meinecke, Phys. Rev. D 91 (2015) 034012, [1410.8063].\n"
         "\n"
         "  For stop-antistop annihilation into gluons and light quarks, please cite:\n"
         "  - M. Klasen, K. Kovařík and L. P. Wiggering, Phys. Rev. D 106 (2022) 115032, [2210.05260].\n"
         "\n");
}


void parseArgs(int argc, char *argv[], Args &args) {
  //Parameters *params = args.params;
  // A value different than "" for an argument means that the user has set
  // that argment through the command-line.
  args.arguments["muR"] = "";
  args.arguments["particleA"] = "";
  args.arguments["particleB"] = "";
  args.arguments["particle1"] = "";
  args.arguments["particle2"] = "";
  args.arguments["output"] = "";
  args.arguments["slha"] = "";
  args.arguments["pcm"] = "";
  args.arguments["renscheme"] = "";
  args.arguments["choosesol"] = "";
  args.arguments["result"] = "";
  args.arguments["formfactor"] = "";
  optind = 1;
  for (;;) {
    static struct option long_options[] = {{"version", no_argument, 0, 'v'},
                                           {"help", no_argument, 0, 'h'},
                                           {"parameter-log", required_argument, 0, 'p'},
                                           {"lo", no_argument, 0, 'l'},
                                           {"nlo", no_argument, 0, 'n'},
                                           {"sommerfeld", no_argument, 0, 's'},
                                           {"full", no_argument, 0, 'f'},
                                           {"result", required_argument, 0, 'r'},
                                           {"muR", required_argument, 0, 'R'},
                                           {"particleA", required_argument, 0, 'a'},
                                           {"particleB", required_argument, 0, 'b'},
                                           {"particle1", required_argument, 0, 'x'},
                                           {"particle2", required_argument, 0, 'y'},
                                           {"output", required_argument, 0, 'o'},
                                           {"slha", required_argument, 0, 'i'},
                                           {"pcm", required_argument, 0, 'e'},
                                           {"renscheme", required_argument, 0, 'z'},
                                           {"choosesol", required_argument, 0, 'c'},
                                           {"legacy", no_argument, 0, 'u'},
                                           {"DD", no_argument, 0, 'd'},
                                           {"formfactor", required_argument, 0, 'F'},
                                           {0, 0, 0, 0}};
    int option_index = 0;
    int c;
    c = getopt_long(argc, argv, "vhp:lnm:t:a:b:f:r:o:", long_options, &option_index);
    if (c == -1) {
      break;
    }
    switch (c) {
    case 'v':
      fprintf(stdout, "DM@NLO " DMNLO_VERSION "\n");
      exit(0);
      break;
    case 'h':
      fprintf(stdout, "Please see https://dmnlo.hepforge.org/ for instructions.");
      exit(0);
      break;
    case 'd':
      args.useDD = true;
      break;
    case 'l':
      args.stopAfterLO = 1;
      break;
    case 'n':
      args.stopAfterNLO = 1;
      break;
    case 's':
      args.SommerfeldAlone = 1;
      break;
    case 'f':
      args.full = 1;
      break;
    case 'p':
      args.logFile = optarg;
      break;
    case 'u':
      args.legacy = 1;
      break;
    case 'r':
      args.arguments["result"] = optarg;
      break;
    case 'R':
      args.arguments["muR"] = optarg;
      break;
    case 'a':
      args.arguments["particleA"] = optarg;
      break;
    case 'b':
      args.arguments["particleB"] = optarg;
      break;
    case 'x':
      args.arguments["particle1"] = optarg;
      break;
    case 'y':
      args.arguments["particle2"] = optarg;
      break;
// LW: For what purpose did I include this option?
    case 'o':
      args.arguments["output"] = optarg;
      break;
    case 'i':
      args.arguments["slha"] = optarg;
      break;
    case 'e':
      args.arguments["pcm"] = optarg;
      break;
    case 'z':
      args.arguments["renscheme"] = optarg;
      break;
    case 'c':
      args.arguments["choosesol"] = optarg;
      break;
    case 'F':
      args.arguments["formfactor"] = optarg;
      break;
    case '?':
      break;
    default:
      break;
    }
  }

  if(args.stopAfterLO + args.stopAfterNLO + args.SommerfeldAlone + args.full > 1) {
    fprintf(stderr, "ERROR: Too many accuracy options were set.\n");
    exit(1);
  } 

  // There can be either 0 or 1 left arguments (input file name).
  if (optind > argc) {
    fprintf(stderr, "ERROR: Too many arguments.\n");
    exit(1);
  } else if (optind == argc - 1) {
    args.inputFile = argv[optind];
  } else {
    fprintf(stderr, "ERROR: No input file was specified.\n");
    exit(1);
  }

  // Read input file
  args.config = readInputFile(args.inputFile);

  // Determine form of output cross section 
  if (get_option("result") == "s") args.fluxFlag = 1;
  else if (get_option("result") == "sv") args.fluxFlag = 0;
  else {
    fprintf(stderr,
            "ERROR: Computation '%s' not recognized. "
            "'sv' or 's' expected.\n",
            get_option("result").c_str());
    exit(1);
  }

  // set PDG numbers of colliding particles
  args.n1 = atol(get_option("particle1").c_str());
  args.n2 = atol(get_option("particle2").c_str());
  args.na = atol(get_option("particleA").c_str());
  args.nb = atol(get_option("particleB").c_str());

  // set renormalization and Coulomb scale, muC is irrelevant for LO and NLO
  args.muR = atof(get_option("muR").c_str());

  // set center-of-mass momentum
  args.pcm = atof(get_option("pcm").c_str());

  // set renormalization scheme
  args.renscheme = atoi(get_option("renscheme").c_str());

  // define solution for mixing matrices
  args.choosesol = atoi(get_option("choosesol").c_str());

  // define the form factor for direct detection
  args.formfactor = atoi(get_option("formfactor").c_str());

  // set SLHA file
  args.slhaFile = get_option("slha").c_str();
}


bool canImprove_ChiChi2QQ(long n1, long n2, long n3, long n4)
{
	//Neu_i + Neu_j -> q qbar
	if((n1 == 1000022 || n1 == 1000023 || n1 == 1000025 || n1 == 1000035) && (n2 == 1000022 || n2 == 1000023 || n2 == 1000025 || n2 == 1000035)){
		if ((n3==1 && n4==-1) || (n3==2 && n4==-2) || (n3==3 && n4==-3) || (n3==4 && n4==-4) || (n3==5 && n4==-5) || (n3==6 && n4==-6) 
			|| (n4==1 && n3==-1) || (n4==2 && n3==-2) || (n4==3 && n3==-3) || (n4==4 && n3==-4) || (n4==5 && n3==-5) || (n4==6 && n3==-6)) {
			return true;			
		}
	} 
	
	//Neu_i + Cha_j -> q q'bar
	if(((n1 == 1000024 || n1 == 1000037 || n1 == -1000024 || n1 == -1000037) && (n2 == 1000022 || n2 == 1000023 || n2 == 1000025 || n2 == 1000035))
	   || ((n2 == 1000024 || n2 == 1000037 || n2 == -1000024 || n2 == -1000037) && (n1 == 1000022 || n1 == 1000023 || n1 == 1000025 || n1 == 1000035))){
		if ((n3==2 && n4==-1) || (n3==4 && n4==-3) || (n3==6 && n4==-5) || (n4==2 && n3==-1) || (n4==4 && n3==-3) || (n4==6 && n3==-5)) {
			return true;
		}
	}
	
	//Cha_i + Cha_j -> q qbar
	if(((n1 == 1000024 || n1 == 1000037) && (n2 == -1000024 || n2 == -1000037)) || ((n1 == -1000024 || n1 == -1000037) && (n2 == 1000024 || n2 == 1000037))){
		if ((n3==1 && n4==-1) || (n3==2 && n4==-2) || (n3==3 && n4==-3) || (n3==4 && n4==-4) || (n3==5 && n4==-5) || (n3==6 && n4==-6) 
			|| (n4==1 && n3==-1) || (n4==2 && n3==-2) || (n4==3 && n3==-3) || (n4==4 && n3==-4) || (n4==5 && n3==-5) || (n4==6 && n3==-6)) {
			return true;
		}
	}
	
	return false;
}

void ChiChi2QQ_setParticleFlags(int *mflags, long n1, long n2, long n3, long n4)
{
	//Neu_i + Neu_j -> q qbar
	if((n1 == 1000022 || n1 == 1000023 || n1 == 1000025 || n1 == 1000035) && (n2 == 1000022 || n2 == 1000023 || n2 == 1000025 || n2 == 1000035)){
		
		mflags[0] = 0;
		
		if(n1 == 1000022) mflags[1] = 1;
		if(n1 == 1000023) mflags[1] = 2;
		if(n1 == 1000025) mflags[1] = 3;
		if(n1 == 1000035) mflags[1] = 4;
		
		if(n2 == 1000022) mflags[2] = 1;
		if(n2 == 1000023) mflags[2] = 2;
		if(n2 == 1000025) mflags[2] = 3;
		if(n2 == 1000035) mflags[2] = 4;
		
		if(n3==1 || n4==1){
			mflags[4] = 1;
			mflags[3] = 4;
		}
		
		if(n3==2 || n4==2){
			mflags[4] = 1;
			mflags[3] = 3;
		}
		
		if(n3==3 || n4==3){
			mflags[4] = 2;
			mflags[3] = 4;
		}
		
		if(n3==4 || n4==4){
			mflags[4] = 2;
			mflags[3] = 3;
		}
		
		if(n3==5 || n4==5){
			mflags[4] = 3;
			mflags[3] = 4;
		}
		
		if(n3==6 || n4==6){
			mflags[4] = 3;
			mflags[3] = 3;
		}
		
	}
	
	//Neu_i + Cha_j -> q q'bar
	if(((n1 == 1000024 || n1 == 1000037 || n1 == -1000024 || n1 == -1000037) && (n2 == 1000022 || n2 == 1000023 || n2 == 1000025 || n2 == 1000035))
	   || ((n2 == 1000024 || n2 == 1000037 || n2 == -1000024 || n2 == -1000037) && (n1 == 1000022 || n1 == 1000023 || n1 == 1000025 || n1 == 1000035))){
		
		mflags[0] = 1;
		mflags[3] = 3;
		
		if(n1 == 1000024 || n1 == -1000024) mflags[1] = 1;
		if(n1 == 1000037 || n1 == -1000037) mflags[1] = 2;
		
		if(n2 == 1000022) mflags[2] = 1;
		if(n2 == 1000023) mflags[2] = 2;
		if(n2 == 1000025) mflags[2] = 3;
		if(n2 == 1000035) mflags[2] = 4;
		
		if(n2 == 1000024 || n2 == -1000024) mflags[1] = 1;
		if(n2 == 1000037 || n2 == -1000037) mflags[1] = 2;
		
		if(n1 == 1000022) mflags[2] = 1;
		if(n1 == 1000023) mflags[2] = 2;
		if(n1 == 1000025) mflags[2] = 3;
		if(n1 == 1000035) mflags[2] = 4;
		
		if(n3==2 || n4 ==2){
			mflags[4] = 1;
		}
		
		if(n3==4 || n4 ==4){
			mflags[4] = 2;
		}
		
		if(n3==6 || n4 ==6){
			mflags[4] = 3;
		}
	}
	
	//Cha_i + Cha_j -> q qbar
	if(((n1 == 1000024 || n1 == 1000037) && (n2 == -1000024 || n2 == -1000037)) || ((n1 == -1000024 || n1 == -1000037) && (n2 == 1000024 || n2 == 1000037))){
		
		mflags[0] = 2;
		
		if(n1 == 1000024) mflags[1] = 1;
		if(n1 == 1000037) mflags[1] = 2;
		
		if(n2 == -1000024) mflags[2] = 1;
		if(n2 == -1000037) mflags[2] = 2;
		
		if(n2 == 1000024) mflags[1] = 1;
		if(n2 == 1000037) mflags[1] = 2;
		
		if(n1 == -1000024) mflags[2] = 1;
		if(n1 == -1000037) mflags[2] = 2;
		
		if(n3==1 || n4==1){
			mflags[4] = 1;
			mflags[3] = 4;
		}
		
		if(n3==2 || n4==2){
			mflags[4] = 1;
			mflags[3] = 3;
		}
		
		if(n3==3 || n4==3){
			mflags[4] = 2;
			mflags[3] = 4;
		}
		
		if(n3==4 || n4==4){
			mflags[4] = 2;
			mflags[3] = 3;
		}
		
		if(n3==5 || n4==5){
			mflags[4] = 3;
			mflags[3] = 4;
		}
		
		if(n3==6 || n4==6){
			mflags[4] = 3;
			mflags[3] = 3;
		}
	}	
}

bool canImprove_NeuQ2qx(long n1, long n2, long n3, long n4)
{
	n1 = labs(n1); n2 = labs(n2); n3 = labs(n3); n4 = labs(n4);
	
	// select all co-annihilation final states
	if(n1==1000022)	// Check if n1 corresponds to a gaugino
	{	if(n2 == 1000006 || n2 == 2000006) 	// Check if n2 corresponds to a Sfermion
		{	if(n3 == 5 || n3 == 6) 										// Check if n3 corresponds to a fermion
			{   if(n4==21 || n4==23 || n4==24 || n4==25 || n4==35 || n4==36 || n4==37)
				{														// Check if n4 corresponds to a gauge or Higgs boson
				 if((n2 > 0 && n3 < 0) || (n2 < 0 && n3 > 0)) return false;
         return true;
				}
			}
			else if(n4 == 5 || n4 == 6) 									// Check if n4 corresponds to a fermion
			{ 	if(n3==21 || n3==23 || n3==24 ||n3==25 || n3==35 || n3==36 || n3==37)
				{														// Check if n3 corresponds to a gauge or Higgs boson
					if((n2 > 0 && n4 < 0) || (n2 < 0 && n4 > 0)) return false;
          return true;
				}
			}	
		}		
	}
	else
	{	if(n2==1000022) 		// Check if n2 corresponds to a gaugino
		{	if(n1 == 1000006 || n1 == 2000006) 	// Check if n1 corresponds to a Sfermion
			{	if(n3 == 5 || n3 == 6) 										// Check if n3 corresponds to a fermion
				{	if (n4==21 || n4==23 || n4==24 || n4==25 || n4==35 || n4==36 || n4==37)
					{														// Check if n4 corresponds to a gauge or Higgs boson
						if((n1 > 0 && n3 < 0) || (n1 < 0 && n3 > 0)) return false;
            return true;
					}
				}
				else if(n4 == 5 || n4 == 6) 									// Check if n4 corresponds to a fermion
				{ 	if (n3 == 21 || n3 == 23 || n3 == 24 ||n3 == 25 || n3 == 35 || n3 == 36 || n3 == 37)
					{														// Check if n4 corresponds to a gauge or Higgs boson
						if((n1 > 0 && n4 < 0) || (n1 < 0 && n4 > 0)) return false;
            return true;
					}
				}	
			}
		}
	}
	return false;
} 

bool canImprove_staustau2QQ(long na, long nb, long n1, long n2){
  if(labs(na) == 1000015 && nb == -na){ // check if na and b correspond to staus
    if(labs(n1) == 6 && n2 == -n1) return true; // check if n1 and n2 correspond to tops
  }
  return false;
}

void staustau2QQ_setParticleFlags(int *mflags, long na, long nb, long n1, long n2){
  if(labs(na) == 1000015 && nb == -na){ // check if na and b correspond to staus
    if(labs(n1) == 6 && n2 == -n1) {
      // indices of na
      mflags[0] = 1;
      mflags[1] = 2;
      mflags[2] = 3;
      // indices of nb
      mflags[3] = 1;
      mflags[4] = 2;
      mflags[5] = 3;
      // indices of n1
      mflags[6] = 3;
      mflags[7] = 3;
      // indices of n2
      mflags[8] = 3;
      mflags[9] = 3;
    }
  }
}


bool canImprove_stst2QQ(long n1, long n2, long n3, long n4){

// Stop_i + Stop_j -> tt    
    if((n1 == 1000006 || n1 == 2000006) && (n2 == 1000006 || n2 == 2000006)){
        if (n3 == 6 && n4 ==6){
            return true;
        }
    }
// Sbottom_i + Sbottom_j -> bb      
    if((n1 == 1000005 || n1 == 2000005) && (n2 == 1000005 || n2 == 2000005)){
        if (n3 == 5 && n4 ==5){
            return true;
        }
    }
// Stop_i + Sbottom_j -> tb  
    if((n1 == 1000006 || n1 == 2000006) && (n2 == 1000005 || n2 == 2000005)){
        if (n3 == 6 && n4 ==5){
            return true;
        }
    }
// Sbottom_i + Stop_j -> bt      
    if((n1 == 1000005 || n1 == 2000005) && (n2 == 1000006 || n2 == 2000006)){
        if (n3 == 5 && n4 ==6){
            return true;
        }
    }
    return false;
}

bool canImprove_stsT2xx(long na, long nb, long n1, long n2){
  if(n1 == 21 && n2 == n1){
    if(abs(na) == 1000006 && nb == -na) return true;
    if(abs(na) == 2000006 && nb == -na) return true;
  } 
  return false;
}
void stst2QQ_setParticleFlags(int *mflags_qq2qq, long n1, long n2, long n3, long n4){
  //Set first flag index of Sfermion1
  if(n1 == 1000005 || n1 == 1000006){
        mflags_qq2qq[0] = 1;
        
  }         
  if(n1 == 2000005 || n1 == 2000006){
        mflags_qq2qq[0] = 2;
        
  }
        
  //Set second flag type if Sfermion1
  if(n1 == 1000006 || n1 == 2000006){
      mflags_qq2qq[1] = 3;
  }
  if(n1 == 1000005 || n1 == 2000005){
      mflags_qq2qq[1] = 4;
  }

  mflags_qq2qq[2] = 3;

  //Set 4 flag index of Sfermion2
  if(n2 == 1000005 || n2 == 1000006){
        mflags_qq2qq[3] = 1;
        
  }         
  if(n2 == 2000005 || n2 == 2000006){
        mflags_qq2qq[3] = 2;
        
  }
        
  //Set second flag type if Sfermion2
  if(n2 == 1000006 || n2 == 2000006){
      mflags_qq2qq[4] = 3;
  }
  if(n2 == 1000005 || n2 == 2000005){
      mflags_qq2qq[4] = 4;
  }

  mflags_qq2qq[5] = 3;

  // set 7th flag, type of fermion1
  if(n3 == 6){
      mflags_qq2qq[6] = 3;
  }
  if(n3 == 5){
      mflags_qq2qq[6] = 4;
  }
  // set 8th flag, gen of fermion1
  mflags_qq2qq[7] = 3;

  // set 9th flag, type of fermion1
  if(n4 == 6){
      mflags_qq2qq[8] = 3;
  }
  if(n4 == 5){
      mflags_qq2qq[8] = 4;
  }
  // set 10th flag, gen of fermion1
  mflags_qq2qq[9] = 3;
}

bool canImprove_QQ2xx(long n1, long n2, long n3, long n4, int ELWDecays)
{	
	// The option ELWdecays = 1 allows to include processes like stop antistop -> gluon + elw. final states
	// These have to be taken into account in the relic routine to avoid double counting, see section III.C of
	// QQ2xx paper (1410.8063). For other purposes e.g. estimating the DMNLO contribution we set ELWdecays = 0.
	
	if(labs(n1) == 1000006 && labs(n2) == 1000006 ){ // stop initial state
		
		if( (labs(n3) >= 22 && labs(n3) <= 37) && (labs(n4) >= 22 && labs(n4) <= 37) ){	  //elw. bosonic final states
			return true;
		}
			
		if( (labs(n3) >= 11 && labs(n3) <= 16) && (labs(n4) >= 11 && labs(n4) <= 16) ){   //elw. fermionic final states
			return true; 
		}
			
		if((ELWDecays==1)&&((labs(n3) == 21 && labs(n4) >= 22 && labs(n4) <= 37) || (labs(n4) == 21 && labs(n3) >= 22 && labs(n3) <= 37))){ //Gluon + elw. final states to correct for double counting in 2->3 processes (on shell decay)
			return true;
		}
		
	}	
	return false;
}


void computeDD(Args &args, double *dd_results){
  ///////////////////////////////   direct detection    //////////////////////////////////
  double pcm, costh, results[8];
  int LoopMode,DenomMode,BoxMode,NuclMode,iflags[3];

  // set kinematics in the nonrelativistic approximation
  pcm = 0.; // Momentum in CMS frame in GeV (pcm != 0 only for testing purposes)	
  costh = 0.5;			// Scattering angle (if pcm = 0, then costheta is irrelevant)	
  
  // LoopMode: 0 -> only tree level, 1 -> Full NLO
  if(args.stopAfterLO) LoopMode = 0;
  else  LoopMode = 1;

  DenomMode = 1; // 0 -> Heavy Squark limit, 1 -> Full kinematics
  BoxMode = 1;  // 0 -> Old result for v != 0, 1 -> New result for v = 0
  NuclMode = 4-args.formfactor; // Nuclear input; 1 -> micrOMEGAs 2.4.1, 2 -> micrOMEGAs 4.0.3, 3 -> DarkSUSY, 4 -> DM@NLO

  iflags[0] = args.renscheme;
  iflags[1] = args.legacy;
  iflags[2] = args.choosesol;

  directdetection_(pcm, costh, LoopMode, DenomMode, BoxMode, NuclMode, args.muR, results, iflags, args.slhaFile.c_str());

  dd_results[0] = results[0]; // g^SI_p
  dd_results[1] = results[2]; // g^SD_p
  dd_results[2] = results[4]; // g^SI_n
  dd_results[3] = results[6]; // g^SD_n
  return;
}


double compute(Args &args, int verbosity){
///////////////////////////////   ChiChi2QQ    //////////////////////////////////
  if(args.correct_ChiChi2QQ && canImprove_ChiChi2QQ(args.na, args.nb, args.n1, args.n2)){
    int iflags[14];
    // use the dipole subtraction method 
    double cut = 0.,result[9];
    // convert PDG numbers to internal convention 
    ChiChi2QQ_setParticleFlags(iflags, args.na, args.nb, args.n1, args.n2);
    //REMAINING CHICHI2QQ DEFAULT SETTINGS
    // diagonalize (0 = Masses from MO, 1 = Diagonalize)
    iflags[5] = 0;
    // iflux (1 = sigma, 0 = v.sigma)
    iflags[6] = args.fluxFlag;
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) iflags[7] = 1;
    else  iflags[7] = 0;
    // Kinematics in MO Cross Section (0 = Ours, 1 = MO)
    iflags[8] = 0;
    // Resummation Style for the Higgs-b-bbar coupling (0 Old DMNLO Style, 1 Spira Style)
    iflags[9] = 1;
    // NNLO Resummation (1 = on, 0 = off), Note: The use NNLO Res. = 1 sets Res. Style to 1 automatically.
    iflags[10] = 1;
    // renormalization scheme
    iflags[11] = args.renscheme;
    // legacy 
    iflags[12]= args.legacy;
    // solution for mixing matrices
    iflags[13] = args.choosesol;

    chichi2qq_(args.muR,args.pcm,iflags,cut,result,args.slhaFile.c_str());
    if(args.stopAfterLO) return result[0];
    else return result[8]; // NLO result, no QED or QCD Sommerfeld enhancement possible
  }
///////////////////////////////   NeuQ2qx    //////////////////////////////////
  else if(args.correct_NeuQ2qx && canImprove_NeuQ2qx(args.na, args.nb, args.n1, args.n2)){
    int iflags[8];
    // use the dipole subtraction method 
    double result[11],error[9];

    // diagonalize (0 = Masses from MO, 1 = Diagonalize)
    iflags[0] = 0;
    // iflux (1=sigma in GeV^-2, 2=sigma in pb, 0=v.sigma)
    iflags[1] = args.fluxFlag;
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) iflags[2] = 1;
    else iflags[2] = 0;
    // Kinematics in MO Cross Section (0 = Ours, 1 = MO)
    iflags[3] = 0;
    // Diagonalized CKM-Matrix = 0, undiagonalized CKM-Matrix = 1
    iflags[4] = 0;
    // renormalization scheme
    iflags[5] = args.renscheme;
    // legacy option (definition of Weinberg angle)
    iflags[6] = args.legacy;
    // solution for mixing matrices
    iflags[7] = args.choosesol;

    neuq2qx_(args.muR, args.na, args.nb, args.n1, args.n2, args.pcm,iflags,result, error, args.slhaFile.c_str());
    if(args.stopAfterLO) return result[0];
    else return result[8]; // NLO result, no QED or QCD Sommerfeld enhancement possible
  }
///////////////////////////////   staustau2QQ    //////////////////////////////////
  else if(args.correct_staustau2QQ && canImprove_staustau2QQ(args.na, args.nb, args.n1, args.n2)){
    int iflags[17];
    double result[9];

    // convert PDG numbers to internal convention 
    staustau2QQ_setParticleFlags(iflags, args.na, args.nb, args.n1, args.n2);
    // diagonalize (0 = Masses from MO, 1 = Diagonalize)
    iflags[10] = 0;
    // iflux (0=v.sigma, 1=sigma in GeV^-2)
    iflags[11] = args.fluxFlag;
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) {iflags[12] = 1; iflags[14] = 0;} 
    else if(args.stopAfterNLO) {iflags[12] = 0; iflags[14] = 0;}
    else if(args.SommerfeldAlone) {iflags[12] = 1; iflags[14] = 1;} 
    else {iflags[12] = 0; iflags[14] = 1;} 
    // renormalization scheme
    iflags[13] = args.renscheme;
    // legacy
    iflags[15] = args.legacy;
    // solution for mixing matrices
    iflags[16] = args.choosesol;

    staustau2qq_(args.muR,args.pcm,iflags,result,args.slhaFile.c_str());

    if(args.stopAfterLO || args.SommerfeldAlone) return result[0];
    else return result[8];
  }
///////////////////////////////   stst2QQ    //////////////////////////////////
  else if(args.correct_stst2QQ && canImprove_stst2QQ(args.na, args.nb, args.n1, args.n2)){
    int iflags[18];
    double result[13],error[13];

    // convert PDG numbers to internal convention 
    stst2QQ_setParticleFlags(iflags, args.na, args.nb, args.n1, args.n2);
    // diagonalize (0 = Masses from MO, 1 = Diagonalize)
    iflags[10] = 0;
    // iflux (0=v.sigma, 1=sigma in GeV^-2)
    iflags[11] = args.fluxFlag;
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) {iflags[12] = 1; iflags[14] = 0;} 
    else if(args.stopAfterNLO) {iflags[12] = 0; iflags[14] = 0;}
    else if(args.SommerfeldAlone) {iflags[12] = 1; iflags[14] = 1;} 
    else {iflags[12] = 0; iflags[14] = 1;} 
    // iQCDEW - 0 - all QCD^2 + QCD-EW + EW^2
    iflags[13] = 0;
    // renormalization scheme
    iflags[15] = args.renscheme;
    // legacy
    iflags[16] = args.legacy;
    // solution for mixing matrices
    iflags[17] = args.choosesol;

    stst2qq_(args.muR,args.pcm,iflags,result,error,args.slhaFile.c_str());

    if(args.stopAfterLO || args.SommerfeldAlone) return result[0];
    else return result[8];
  }
///////////////////////////////   QQ2xx    //////////////////////////////////
  else if(args.correct_QQ2xx && canImprove_QQ2xx(args.na, args.nb, args.n1, args.n2, 1)){
    int iflags[8];
    double result[10],error[8];
    long n5 = 21;

    iflags[0] = 0;                   // diagonalize
    iflags[1] = args.fluxFlag;       // iflux (1=sigma, 0=v.sigma)
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) {iflags[2] = 1; iflags[5] = 0;} 
    else if(args.stopAfterNLO) {iflags[2] = 0; iflags[5] = 0;}
    else if(args.SommerfeldAlone) {iflags[2] = 1; iflags[5] = 1;} 
    else {iflags[2] = 0; iflags[5] = 1;} 
    iflags[3] = 0;       // kinematics in MO (0=our, 1=MO original)  
    iflags[4] = args.renscheme;   // renormalization scheme
    // legacy
    iflags[6] = args.legacy;
    // choosesol
    iflags[7] = args.choosesol;

    qq2xx_(args.muR,args.na, args.nb, args.n1, args.n2,n5,args.pcm,iflags,result,error,args.slhaFile.c_str());

    if(args.stopAfterLO) return result[0];
    else if(args.SommerfeldAlone) return result[6];
    else return result[1];
  }
///////////////////////////////   stsT2xx    //////////////////////////////////
  else if(args.correct_stsT2xx && canImprove_stsT2xx(args.na, args.nb, args.n1, args.n2)){
    int iflags[13] = {1,3,3,1,3,3,0,args.fluxFlag,1,0,args.renscheme,args.legacy,args.choosesol};
    double result[15], error[15], gaugek1 = 0., gaugek2 = 0.; 
    // only tree (1 = only tree, 0 = Loops)
    if(args.stopAfterLO) {iflags[8] = 1; iflags[9] = 0;} 
    else if(args.stopAfterNLO) {iflags[8] = 0; iflags[9] = 0;}
    else if(args.SommerfeldAlone) {iflags[8] = 1; iflags[9] = 1;} 
    else {iflags[8] = 0; iflags[9] = 1;} 

    stst2xx_(args.muR,args.pcm,iflags,result,error,gaugek1,gaugek2,args.slhaFile.c_str());

    if(args.stopAfterLO) {
      if(verbosity==0) printf("~t1 ~t1 -> gg: %e GeV^-2\n~t1 ~t1 -> qqbar: %e GeV^-2\n",result[0],result[1]);
      return result[0]+result[1];
    }
    else if(args.SommerfeldAlone) return result[13];
    else return result[0]+result[1]+result[2]+result[6]+result[10]+result[11]+result[12]+result[13];

  }
  // set cross section with light quarks in the final state to zero since this contribution is already accounted for in stsT2xx.
  // the return value is with 1e-200 slightly larger than zero to ensure that the DM@NLO result for light quarks
  // actually replaces the CalcHEP cross section in the MicrOMEGAs function improveCrossSection.
  else if(args.correct_stsT2xx && (args.n1 == 1 || args.n1 == 2 || args.n1 == 3 || args.n1 == 4) && args.n1 == -args.n2){
    if(args.na == 1000006 || args.na == 2000006) return 1e-200;
  }
  if(verbosity==0) {printf("No correction available.\n"); exit(0);}
	return 0.;

}

int canImprove_dmnlo(long n1, long n2, long n3, long n4){
  if(canImprove_ChiChi2QQ(n1, n2, n3, n4) ||  canImprove_NeuQ2qx(n1, n2, n3, n4) ||
      canImprove_staustau2QQ(n1, n2, n3, n4) || canImprove_stst2QQ(n1, n2, n3, n4) ||
      canImprove_QQ2xx(n1, n2, n3, n4, 0) || canImprove_stsT2xx(n1, n2, n3, n4)) return 1;
	return 0;
}

int consistent_RS_dmnlo(int renscheme, char *slha, double muR){
  return check_renscheme_(renscheme, muR, slha);
}

// Returns the cross section in GeV^-2 and not the cross section times velocity
double cs_dmnlo(int order,long na, long nb, long n1, long n2, double pcm, double muR, char *slha, int renscheme, int choosesol, int *corrFlags){
  struct Args args;
  double result;
  args.na = na; args.nb = nb; args.n1 = n1; args.n2 = n2;
  args.pcm = pcm; args.muR = muR; args.slhaFile = slha; args.renscheme = renscheme; args.choosesol = choosesol;
  args.fluxFlag = 1;

// --- specify the order/accuracy of the calculation
  if(order==0) args.stopAfterLO = 1; // LO
  if(order==1)  args.stopAfterNLO = 1; // NLO
  if(order==2)  {args.stopAfterLO = 0; args.stopAfterNLO = 0;} // Full
  if(order==3)  args.SommerfeldAlone = 1; // Sommerfeld alone

// --- decide which processes are included
  if(corrFlags[0] == 0) args.correct_ChiChi2QQ = false;
  if(corrFlags[1] == 0) args.correct_NeuQ2qx = false;
  if(corrFlags[2] == 0) args.correct_staustau2QQ = false;
  if(corrFlags[3] == 0) args.correct_stst2QQ= false;
  if(corrFlags[4] == 0) args.correct_QQ2xx = false;
  if(corrFlags[5] == 0) args.correct_stsT2xx = false;

  result = compute(args,1);

  return result;
}

// function for the direct detection module
void dd_dmnlo(int order, double muR, char *slha, int renscheme, int choosesol, int formfactor, double *cs){
  struct Args args;
  args.muR = muR; args.slhaFile = slha; args.renscheme = renscheme; args.choosesol = choosesol; args.formfactor = formfactor;

// --- specify the order/accuracy of the calculation (LO vs. NLO)
  if(order==0) args.stopAfterLO = 1; // LO
  if(order==1)  args.stopAfterNLO = 1; // NLO

// --- compute the SI and SD DM-proton and DM-neutron cross sections
// cs[0] = SI proton; cs[1] = SD proton; cs[2] = SI neutron; cs[3] = SD neutron; 
  computeDD(args,cs);

  return;
}