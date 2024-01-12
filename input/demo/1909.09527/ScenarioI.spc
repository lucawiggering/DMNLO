# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 14.08.2023,  12:05
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.3.3      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
# Either the general MSSM or a model has been used
# which has not yet been implemented in the LesHouches standard
Block MINPAR  # Input parameters
    3    1.88993218E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    5.45980000E+02  # M_1
    2   -3.46166000E+03  # M_2
    3    3.12637000E+03  # M_3
   11   -5.24671000E+03  # A_t
   12   -2.53026000E+03  # A_b
   13    1.58636000E+03  # A_tau
   23    2.62335000E+03  # mu
   25    1.80400000E+01  # tan(beta)
   26    2.96227000E+03  # m_A, pole mass
   31    3.58673000E+03  # M^2_L11
   32    3.58673000E+03  # M^2_L22
   33    5.62570000E+02  # M^2_L33
   34    3.95037000E+03  # M^2_E11
   35    3.95037000E+03  # M^2_E22
   36    5.85450000E+02  # M^2_E33
   41    1.59946000E+03  # M^2_Q11
   42    1.59946000E+03  # M^2_Q22
   43    3.00698000E+03  # M^2_Q33
   44    3.90439000E+03  # M^2_U11
   45    3.90439000E+03  # M^2_U22
   46    3.09304000E+03  # M^2_U33
   47    3.09670000E+03  # M^2_D11
   48    3.09670000E+03  # M^2_D22
   49    5.81560000E+02  # M^2_D33
Block SMINPUTS  # SM parameters
         1     1.27908953E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637000E-05  # G_mu [GeV^-2]
         3     1.18400000E-01  # alpha_s(MZ)^MSbar
         4     9.11876000E+01  # m_Z(pole)
         5     4.18000000E+00  # m_b(m_b), MSbar
         6     1.73200000E+02  # m_t(pole)
         7     1.77700000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998930E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658372E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     2.50000000E-03  # m_u(2 GeV), MSbar
        23     9.50000000E-02  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  3.04733379E+03  # (SUSY scale)
   1    3.64616115E-01  # g'(Q)^DRbar
   2    6.34618746E-01  # g(Q)^DRbar
   3    1.00936520E+00  # g3(Q)^DRbar
Block Yu Q=  3.04733379E+03  # (SUSY scale)
  1  1     6.80932905E-06   # Y_u(Q)^DRbar
  2  2     3.45929002E-03   # Y_c(Q)^DRbar
  3  3     8.22588681E-01   # Y_t(Q)^DRbar
Block Yd Q=  3.04733379E+03  # (SUSY scale)
  1  1     2.00830126E-04   # Y_d(Q)^DRbar
  2  2     3.81577049E-03   # Y_s(Q)^DRbar
  3  3     2.01277300E-01   # Y_b(Q)^DRbar
Block Ye Q=  3.04733379E+03  # (SUSY scale)
  1  1     5.11973266E-05   # Y_e(Q)^DRbar
  2  2     1.05859770E-02   # Y_mu(Q)^DRbar
  3  3     1.70308653E-01   # Y_tau(Q)^DRbar
Block Au Q=  3.04733379E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3    -5.24671000E+03   # A_t(Q)^DRbar
Block Ad Q=  3.04733379E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3    -2.53026000E+03   # A_b(Q)^DRbar
Block Ae Q=  3.04733379E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     1.58636000E+03   # A_tau(Q)^DRbar
Block MSOFT Q=  3.04733379E+03  # soft SUSY breaking masses at Q
   1    5.45980000E+02  # M_1
   2   -3.46166000E+03  # M_2
   3    3.12637000E+03  # M_3
  21    1.72938159E+06  # M^2_(H,d)
  22   -6.69629679E+06  # M^2_(H,u)
  31    3.58673000E+03  # M_(L,11)
  32    3.58673000E+03  # M_(L,22)
  33    5.62570000E+02  # M_(L,33)
  34    3.95037000E+03  # M_(E,11)
  35    3.95037000E+03  # M_(E,22)
  36    5.85450000E+02  # M_(E,33)
  41    1.59946000E+03  # M_(Q,11)
  42    1.59946000E+03  # M_(Q,22)
  43    3.00698000E+03  # M_(Q,33)
  44    3.90439000E+03  # M_(U,11)
  45    3.90439000E+03  # M_(U,22)
  46    3.09304000E+03  # M_(U,33)
  47    3.09670000E+03  # M_(D,11)
  48    3.09670000E+03  # M_(D,22)
  49    5.81560000E+02  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73200000E+02  # m_t(pole)
        23     9.11876000E+01  # m_Z(pole)
        24     8.04084165E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.28311554E+02  # h0
        35     2.96225612E+03  # H0
        36     2.96227000E+03  # A0
        37     2.96369622E+03  # H+
   1000001     1.79656145E+03  # ~d_L
   2000001     3.21166917E+03  # ~d_R
   1000002     1.79482870E+03  # ~u_L
   2000002     3.99093231E+03  # ~u_R
   1000003     1.79656118E+03  # ~s_L
   2000003     3.21166828E+03  # ~s_R
   1000004     1.79482867E+03  # ~c_L
   2000004     3.99093182E+03  # ~c_R
   1000005     8.50201290E+02  # ~b_1
   2000005     3.10442109E+03  # ~b_2
   1000006     2.98580110E+03  # ~t_1
   2000006     3.22357354E+03  # ~t_2
   1000011     3.60865081E+03  # ~e_L-
   2000011     3.95608693E+03  # ~e_R-
   1000012     3.60737576E+03  # ~nu_eL
   1000013     3.60864652E+03  # ~mu_L-
   2000013     3.95608049E+03  # ~mu_R-
   1000014     3.60737259E+03  # ~nu_muL
   1000015     5.42390534E+02  # ~tau_1-
   2000015     6.67344695E+02  # ~tau_2-
   1000016     6.13880125E+02  # ~nu_tauL
   1000021     3.26944959E+03  # ~g
   1000022     5.40555256E+02  # ~chi_10
   1000023    -2.63259065E+03  # ~chi_20
   1000025     2.63623240E+03  # ~chi_30
   1000035    -3.42900879E+03  # ~chi_40
   1000024     2.63304387E+03  # ~chi_1+
   1000037     3.42908515E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -5.54892292E-02   # alpha
Block Hmix Q=  3.04733379E+03  # Higgs mixing parameters
   1    2.62335000E+03  # mu
   2    1.80400000E+01  # tan[beta](Q)
   3    2.43607900E+02  # v(Q)
   4    8.71834231E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     7.17832882E-01   # Re[R_st(1,1)]
   1  2     6.96215451E-01   # Re[R_st(1,2)]
   2  1    -6.96215451E-01   # Re[R_st(2,1)]
   2  2     7.17832882E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     1.09373547E-02   # Re[R_sb(1,1)]
   1  2     9.99940185E-01   # Re[R_sb(1,2)]
   2  1    -9.99940185E-01   # Re[R_sb(2,1)]
   2  2     1.09373547E-02   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.46821003E-01   # Re[R_sta(1,1)]
   1  2     7.62641849E-01   # Re[R_sta(1,2)]
   2  1    -7.62641849E-01   # Re[R_sta(2,1)]
   2  2     6.46821003E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.99828532E-01   # Re[N(1,1)]
   1  2    -1.29774931E-04   # Re[N(1,2)]
   1  3    -1.79123187E-02   # Re[N(1,3)]
   1  4     4.69446923E-03   # Re[N(1,4)]
   2  1     9.32992847E-03   # Re[N(2,1)]
   2  2     6.41066639E-02   # Re[N(2,2)]
   2  3    -7.06057110E-01   # Re[N(2,3)]
   2  4    -7.05185539E-01   # Re[N(2,4)]
   3  1     1.59851750E-02   # Re[N(3,1)]
   3  2    -9.52827473E-03   # Re[N(3,2)]
   3  3    -7.06874843E-01   # Re[N(3,3)]
   3  4     7.07093801E-01   # Re[N(3,4)]
   4  1    -5.76764745E-04   # Re[N(4,1)]
   4  2     9.97897555E-01   # Re[N(4,2)]
   4  3     3.86065116E-02   # Re[N(4,3)]
   4  4     5.20545273E-02   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -5.49289430E-02   # Re[U(1,1)]
   1  2     9.98490266E-01   # Re[U(1,2)]
   2  1    -9.98490266E-01   # Re[U(2,1)]
   2  2    -5.49289430E-02   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1     7.40277254E-02   # Re[V(1,1)]
   1  2     9.97256184E-01   # Re[V(1,2)]
   2  1     9.97256184E-01   # Re[V(2,1)]
   2  2    -7.40277254E-02   # Re[V(2,2)]
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   1000011     5.00232715E+00   # ~e^-_L
#    BR                NDA      ID1      ID2
     9.11759348E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     3.05192266E-03    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     2.71066000E-02    2     1000035        11   # BR(~e^-_L -> chi^0_4 e^-)
     3.81334909E-03    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
     5.42687019E-02    2    -1000037        12   # BR(~e^-_L -> chi^-_2 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000011     2.01477503E+01   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99889874E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000013     5.00414218E+00   # ~mu^-_L
#    BR                NDA      ID1      ID2
     9.11439677E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     3.22857159E-03    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     1.74045515E-04    2     1000025        13   # BR(~mu^-_L -> chi^0_3 mu^-)
     2.70953465E-02    2     1000035        13   # BR(~mu^-_L -> chi^0_4 mu^-)
     3.81616828E-03    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     5.42461913E-02    2    -1000037        14   # BR(~mu^-_L -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     2.01542396E+01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99563098E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     1.49873376E-04    2     1000025        13   # BR(~mu^-_R -> chi^0_3 mu^-)
     1.33735512E-04    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
DECAY   1000015     3.79880454E-05   # ~tau^-_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
DECAY   2000015     5.19078422E-01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     4.45942157E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     5.53136008E-01    2     1000015        23   # BR(~tau^-_2 -> ~tau^-_1 Z)
#    BR                NDA      ID1      ID2       ID3
     4.61420267E-04    3     1000016         1        -2   # BR(~tau^-_2 -> ~nu_tau d u_bar)
     4.60414405E-04    3     1000016         3        -4   # BR(~tau^-_2 -> ~nu_tau s c_bar)
DECAY   1000012     5.00271793E+00   # ~nu_e
#    BR                NDA      ID1      ID2
     9.10513258E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     2.17815016E-03    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     2.19606526E-04    2     1000025        12   # BR(~nu_e -> chi^0_3 nu_e)
     2.67752825E-02    2     1000035        12   # BR(~nu_e -> chi^0_4 nu_e)
     6.91204442E-03    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
     5.34016581E-02    2     1000037        11   # BR(~nu_e -> chi^+_2 e^-)
DECAY   1000014     5.00445149E+00   # ~nu_mu
#    BR                NDA      ID1      ID2
     9.10196981E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     2.17738502E-03    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     2.19529377E-04    2     1000025        14   # BR(~nu_mu -> chi^0_3 nu_mu)
     2.67651034E-02    2     1000035        14   # BR(~nu_mu -> chi^0_4 nu_mu)
     7.25966730E-03    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     5.33813343E-02    2     1000037        13   # BR(~nu_mu -> chi^+_2 mu^-)
DECAY   1000016     4.30373640E-02   # ~nu_tau
#    BR                NDA      ID1      ID2
     9.50981657E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
     2.45051835E-02    3     1000015        -2         1   # BR(~nu_tau -> ~tau^-_1 u_bar d)
     2.45131594E-02    3     1000015        -4         3   # BR(~nu_tau -> ~tau^-_1 c_bar s)
DECAY   1000001     2.17976497E-01   # ~d_L
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         1   # BR(~d_L -> chi^0_1 d)
DECAY   2000001     1.78167282E+00   # ~d_R
#    BR                NDA      ID1      ID2
     9.99959636E-01    2     1000022         1   # BR(~d_R -> chi^0_1 d)
DECAY   1000003     2.17976870E-01   # ~s_L
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         3   # BR(~s_L -> chi^0_1 s)
DECAY   2000003     1.78271257E+00   # ~s_R
#    BR                NDA      ID1      ID2
     9.99375918E-01    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     2.34406881E-04    2     1000004       -24   # BR(~s_R -> ~c_L W^-)
     1.20713484E-04    2     1000003        23   # BR(~s_R -> ~s_L Z)
     1.18904037E-04    2     1000003        25   # BR(~s_R -> ~s_L h^0)
DECAY   1000005     1.77361417E-01   # ~b_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
DECAY   2000005     6.97283957E+00   # ~b_2
#    BR                NDA      ID1      ID2
     6.13669381E-02    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     1.46432121E-02    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     1.39958933E-02    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     4.20433597E-01    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     2.29201208E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
     1.29724012E-01    2     1000005        23   # BR(~b_2 -> ~b_1 Z)
     1.30635139E-01    2     1000005        25   # BR(~b_2 -> ~b_1 h^0)
DECAY   1000002     2.18273450E-01   # ~u_L
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         2   # BR(~u_L -> chi^0_1 u)
DECAY   2000002     3.23709633E+01   # ~u_R
#    BR                NDA      ID1      ID2
     2.79213482E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     7.20754929E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000004     2.18273241E-01   # ~c_L
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~c_L -> chi^0_1 c)
DECAY   2000004     3.23714991E+01   # ~c_R
#    BR                NDA      ID1      ID2
     2.79208763E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     7.20740811E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000006     6.81793401E+00   # ~t_1
#    BR                NDA      ID1      ID2
     5.07399285E-01    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     1.74497881E-01    2     1000023         6   # BR(~t_1 -> chi^0_2 t)
     6.37687875E-02    2     1000025         6   # BR(~t_1 -> chi^0_3 t)
     1.32588888E-01    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
     1.15383730E-01    2     1000005        24   # BR(~t_1 -> ~b_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     2.55701813E-04    3     1000016       -15         5   # BR(~t_1 -> ~nu_tau tau^+ b)
     3.38046811E-04    3    -1000015        16         5   # BR(~t_1 -> ~tau^+_1 nu_tau b)
     4.31895633E-04    3    -2000015        16         5   # BR(~t_1 -> ~tau^+_2 nu_tau b)
     5.33570279E-03    3     1000022        24         5   # BR(~t_1 -> chi^0_1 W^+ b)
DECAY   2000006     1.98932336E+01   # ~t_2
#    BR                NDA      ID1      ID2
     1.90383751E-01    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     8.77173800E-02    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     1.45664254E-01    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     1.45346731E-01    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     4.86186994E-02    2     1000005        24   # BR(~t_2 -> ~b_1 W^+)
     7.76164089E-02    2     2000005        24   # BR(~t_2 -> ~b_2 W^+)
     3.04351898E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     3.00877866E-04    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     5.95066273E+00   # chi^+_1
#    BR                NDA      ID1      ID2
     9.35330449E-02    2    -1000015        16   # BR(chi^+_1 -> ~tau^+_1 nu_tau)
     2.68126905E-02    2    -2000015        16   # BR(chi^+_1 -> ~tau^+_2 nu_tau)
     1.22414559E-01    2     1000016       -15   # BR(chi^+_1 -> ~nu_tau tau^+)
     8.35228938E-03    2     1000002        -1   # BR(chi^+_1 -> ~u_L d_bar)
     8.40687325E-03    2     1000004        -3   # BR(chi^+_1 -> ~c_L s_bar)
     4.58304207E-03    2    -1000001         2   # BR(chi^+_1 -> ~d^*_L u)
     4.62981466E-03    2    -1000003         4   # BR(chi^+_1 -> ~s^*_L c)
     4.35426801E-01    2    -1000005         6   # BR(chi^+_1 -> ~b^*_1 t)
     2.85106215E-01    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     1.07306445E-02    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
DECAY   1000037     1.29710093E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     4.05439737E-02    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     5.82815758E-02    2    -2000015        16   # BR(chi^+_2 -> ~tau^+_2 nu_tau)
     9.87096325E-02    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     1.66565296E-01    2     1000002        -1   # BR(chi^+_2 -> ~u_L d_bar)
     1.66565328E-01    2     1000004        -3   # BR(chi^+_2 -> ~c_L s_bar)
     1.13885915E-02    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     1.68184359E-03    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     1.66734453E-01    2    -1000001         2   # BR(chi^+_2 -> ~d^*_L u)
     1.66734318E-01    2    -1000003         4   # BR(chi^+_2 -> ~s^*_L c)
     7.98492615E-03    2    -2000005         6   # BR(chi^+_2 -> ~b^*_2 t)
     1.24051757E-04    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     2.71848952E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.75729395E-02    2     1000025        24   # BR(chi^+_2 -> chi^0_3 W^+)
     2.85780773E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     2.99156437E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     4.68226313E-04    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     4.61481207E-04    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
     4.54669038E-04    3     1000025        -5         6   # BR(chi^+_2 -> chi^0_3 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     5.95805801E+00   # chi^0_2
#    BR                NDA      ID1      ID2
     7.24735492E-02    2     1000015       -15   # BR(chi^0_2 -> ~tau^-_1 tau^+)
     7.24735492E-02    2    -1000015        15   # BR(chi^0_2 -> ~tau^+_1 tau^-)
     4.58805249E-02    2     2000015       -15   # BR(chi^0_2 -> ~tau^-_2 tau^+)
     4.58805249E-02    2    -2000015        15   # BR(chi^0_2 -> ~tau^+_2 tau^-)
     2.73131673E-03    2     1000016       -16   # BR(chi^0_2 -> ~nu_tau nu_bar_tau)
     2.73131673E-03    2    -1000016        16   # BR(chi^0_2 -> ~nu^*_tau nu_tau)
     3.30209799E-03    2     1000002        -2   # BR(chi^0_2 -> ~u_L u_bar)
     3.30209799E-03    2    -1000002         2   # BR(chi^0_2 -> ~u^*_L u)
     3.32554617E-03    2     1000004        -4   # BR(chi^0_2 -> ~c_L c_bar)
     3.32554617E-03    2    -1000004         4   # BR(chi^0_2 -> ~c^*_L c)
     2.94381488E-03    2     1000001        -1   # BR(chi^0_2 -> ~d_L d_bar)
     2.94381488E-03    2    -1000001         1   # BR(chi^0_2 -> ~d^*_L d)
     2.97111878E-03    2     1000003        -3   # BR(chi^0_2 -> ~s_L s_bar)
     2.97111878E-03    2    -1000003         3   # BR(chi^0_2 -> ~s^*_L s)
     2.14542133E-01    2     1000005        -5   # BR(chi^0_2 -> ~b_1 b_bar)
     2.14542133E-01    2    -1000005         5   # BR(chi^0_2 -> ~b^*_1 b)
     2.19212040E-01    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     7.74348201E-02    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     6.98646668E-03    3     1000022         6        -6   # BR(chi^0_2 -> chi^0_1 t t_bar)
DECAY   1000025     5.75493508E+00   # chi^0_3
#    BR                NDA      ID1      ID2
     5.66410836E-02    2     1000015       -15   # BR(chi^0_3 -> ~tau^-_1 tau^+)
     5.66410836E-02    2    -1000015        15   # BR(chi^0_3 -> ~tau^+_1 tau^-)
     6.18537273E-02    2     2000015       -15   # BR(chi^0_3 -> ~tau^-_2 tau^+)
     6.18537273E-02    2    -2000015        15   # BR(chi^0_3 -> ~tau^+_2 tau^-)
     2.87393106E-04    2     1000016       -16   # BR(chi^0_3 -> ~nu_tau nu_bar_tau)
     2.87393106E-04    2    -1000016        16   # BR(chi^0_3 -> ~nu^*_tau nu_tau)
     1.25227906E-04    2     1000001        -1   # BR(chi^0_3 -> ~d_L d_bar)
     1.25227906E-04    2    -1000001         1   # BR(chi^0_3 -> ~d^*_L d)
     1.53615275E-04    2     1000003        -3   # BR(chi^0_3 -> ~s_L s_bar)
     1.53615275E-04    2    -1000003         3   # BR(chi^0_3 -> ~s^*_L s)
     2.21932322E-01    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     2.21932322E-01    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     7.83896925E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     2.32678209E-01    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     6.73860502E-03    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
DECAY   1000035     1.29809348E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     2.02268285E-02    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     2.02268285E-02    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     2.90677505E-02    2     2000015       -15   # BR(chi^0_4 -> ~tau^-_2 tau^+)
     2.90677505E-02    2    -2000015        15   # BR(chi^0_4 -> ~tau^+_2 tau^-)
     4.93997334E-02    2     1000016       -16   # BR(chi^0_4 -> ~nu_tau nu_bar_tau)
     4.93997334E-02    2    -1000016        16   # BR(chi^0_4 -> ~nu^*_tau nu_tau)
     8.33029442E-02    2     1000002        -2   # BR(chi^0_4 -> ~u_L u_bar)
     8.33029442E-02    2    -1000002         2   # BR(chi^0_4 -> ~u^*_L u)
     8.33028803E-02    2     1000004        -4   # BR(chi^0_4 -> ~c_L c_bar)
     8.33028803E-02    2    -1000004         4   # BR(chi^0_4 -> ~c^*_L c)
     4.94309655E-03    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     4.94309655E-03    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     3.84154972E-04    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     3.84154972E-04    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     8.32183761E-02    2     1000001        -1   # BR(chi^0_4 -> ~d_L d_bar)
     8.32183761E-02    2    -1000001         1   # BR(chi^0_4 -> ~d^*_L d)
     8.32183889E-02    2     1000003        -3   # BR(chi^0_4 -> ~s_L s_bar)
     8.32183889E-02    2    -1000003         3   # BR(chi^0_4 -> ~s^*_L s)
     5.14795520E-03    2     2000005        -5   # BR(chi^0_4 -> ~b_2 b_bar)
     5.14795520E-03    2    -2000005         5   # BR(chi^0_4 -> ~b^*_2 b)
     2.80604917E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     2.80604917E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     6.26823326E-04    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     2.74185924E-02    2     1000025        23   # BR(chi^0_4 -> chi^0_3 Z)
     2.92396579E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
     6.13369315E-04    2     1000025        25   # BR(chi^0_4 -> chi^0_3 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.86277653E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     2.77017637E-04    3     1000025         6        -6   # BR(chi^0_4 -> chi^0_3 t t_bar)
     4.60569245E-04    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     4.60569245E-04    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     1.88253994E+02   # ~g
#    BR                NDA      ID1      ID2
     8.59062097E-02    2     1000002        -2   # BR(~g -> ~u_L u_bar)
     8.59062097E-02    2    -1000002         2   # BR(~g -> ~u^*_L u)
     8.59061937E-02    2     1000004        -4   # BR(~g -> ~c_L c_bar)
     8.59061937E-02    2    -1000004         4   # BR(~g -> ~c^*_L c)
     1.45453673E-03    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     1.45453673E-03    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     8.57630984E-02    2     1000001        -1   # BR(~g -> ~d_L d_bar)
     8.57630984E-02    2    -1000001         1   # BR(~g -> ~d^*_L d)
     2.16018065E-04    2     2000001        -1   # BR(~g -> ~d_R d_bar)
     2.16018065E-04    2    -2000001         1   # BR(~g -> ~d^*_R d)
     8.57631167E-02    2     1000003        -3   # BR(~g -> ~s_L s_bar)
     8.57631167E-02    2    -1000003         3   # BR(~g -> ~s^*_L s)
     2.16024533E-04    2     2000003        -3   # BR(~g -> ~s_R s_bar)
     2.16024533E-04    2    -2000003         3   # BR(~g -> ~s^*_R s)
     1.52997460E-01    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     1.52997460E-01    2    -1000005         5   # BR(~g -> ~b^*_1 b)
     1.70477716E-03    2     2000005        -5   # BR(~g -> ~b_2 b_bar)
     1.70477716E-03    2    -2000005         5   # BR(~g -> ~b^*_2 b)
#    BR                NDA      ID1      ID2       ID3
DECAY        25     3.14015628E-03   # h^0
#    BR                NDA      ID1      ID2
     2.96389007E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     7.63516801E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     1.60576072E-04    2           3        -3   # BR(h^0 -> s s_bar)
     4.30796856E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.95609184E-02    2           4        -4   # BR(h^0 -> c c_bar)
     1.06135794E-01    2          21        21   # BR(h^0 -> g g)
     3.91557977E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.52858006E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.52858006E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.52858006E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     5.35003020E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     5.35003020E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.52858006E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.52858006E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.52858006E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     5.35003020E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     5.35003020E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     1.24540382E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     1.24540382E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     1.25281694E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     7.41311800E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     5.78223204E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     5.78223204E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     5.59690409E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     4.29960844E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     4.44787080E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     1.09526003E+01   # H^0
#    BR                NDA      ID1      ID2
     6.01117676E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     1.55585851E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     2.69117863E-04    2           3        -3   # BR(H^0 -> s s_bar)
     7.48805178E-01    2           5        -5   # BR(H^0 -> b b_bar)
     3.82610087E-02    2           6        -6   # BR(H^0 -> t t_bar)
     2.37536395E-02    2     1000015  -1000015   # BR(H^0 -> ~tau^-_1 ~tau^+_1)
     6.45905872E-04    2     1000015  -2000015   # BR(H^0 -> ~tau^-_1 ~tau^+_2)
     6.45905872E-04    2    -1000015   2000015   # BR(H^0 -> ~tau^+_1 ~tau^-_2)
     2.34371601E-02    2     2000015  -2000015   # BR(H^0 -> ~tau^-_2 ~tau^+_2)
     1.79172244E-04    2     1000022   1000022   # BR(H^0 -> chi^0_1 chi^0_1)
     2.67361199E-03    2          23        23   # BR(H^0 -> Z Z)
     5.00548262E-03    2          24       -24   # BR(H^0 -> W^+ W^-)
DECAY        36     1.08692774E+01   # A^0
#    BR                NDA      ID1      ID2
     6.05736263E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.56781494E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     2.71185599E-04    2           3        -3   # BR(A^0 -> s s_bar)
     7.54557743E-01    2           5        -5   # BR(A^0 -> b b_bar)
     3.86703861E-02    2           6        -6   # BR(A^0 -> t t_bar)
     2.44340195E-02    2     1000015  -2000015   # BR(A^0 -> ~tau^-_1 ~tau^+_2)
     2.44340195E-02    2    -1000015   2000015   # BR(A^0 -> ~tau^+_1 ~tau^-_2)
     2.20770047E-04    2     1000022   1000022   # BR(A^0 -> chi^0_1 chi^0_1)
DECAY        37     9.70426180E+00   # H^+
#    BR                NDA      ID1      ID2
     6.78782733E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.75687964E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.65246632E-04    2          -3         4   # BR(H^+ -> s_bar c)
     7.68619140E-01    2          -5         6   # BR(H^+ -> b_bar t)
     3.14518216E-02    2    -1000015   1000016   # BR(H^+ -> ~tau^+_1 ~nu_tau)
     2.32856679E-02    2    -2000015   1000016   # BR(H^+ -> ~tau^+_2 ~nu_tau)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00410030E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    3.25437500E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    3.25441600E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99987401E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    3.08534711E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    3.07274792E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00410030E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    3.25437500E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    3.25441600E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999987E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    1.28887937E-08        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999987E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    1.28887937E-08        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15565942E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    1.30739261E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    9.92132306E-02        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    1.28887937E-08        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999987E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.20579191E-04   # BR(b -> s gamma)
    2    1.59041935E-06   # BR(b -> s mu+ mu-)
    3    3.52784222E-05   # BR(b -> s nu nu)
    4    1.94169348E-15   # BR(Bd -> e+ e-)
    5    8.29469042E-11   # BR(Bd -> mu+ mu-)
    6    1.73119870E-08   # BR(Bd -> tau+ tau-)
    7    6.52043434E-14   # BR(Bs -> e+ e-)
    8    2.78552633E-09   # BR(Bs -> mu+ mu-)
    9    5.89172788E-07   # BR(Bs -> tau+ tau-)
   10    9.67291987E-05   # BR(B_u -> tau nu)
   11    9.99190496E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41380007E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93438046E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15605623E-03   # epsilon_K
   17    2.28161062E-15   # Delta(M_K)
   18    2.48212256E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.29457080E-11   # BR(K^+ -> pi^+ nu nu)
   20   -4.09120766E-16   # Delta(g-2)_electron/2
   21   -1.74911582E-11   # Delta(g-2)_muon/2
   22    3.32312006E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    3.63711385E-05   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88850708E-01   # C7
     0305 4422   00   2    -1.88671117E-01   # C7
     0305 4322   00   2     3.34420294E-05   # C7'
     0305 6421   00   0    -9.52378784E-02   # C8
     0305 6421   00   2    -9.84204321E-02   # C8
     0305 6321   00   2    -3.50266536E-05   # C8'
 03051111 4133   00   0     1.16157211E+00   # C9 e+e-
 03051111 4133   00   2     1.16177370E+00   # C9 e+e-
 03051111 4233   00   2     3.10740051E-06   # C9' e+e-
 03051111 4137   00   0    -3.98586015E+00   # C10 e+e-
 03051111 4137   00   2    -3.98476007E+00   # C10 e+e-
 03051111 4237   00   2    -8.83778658E-05   # C10' e+e-
 03051313 4133   00   0     1.16157211E+00   # C9 mu+mu-
 03051313 4133   00   2     1.16177360E+00   # C9 mu+mu-
 03051313 4233   00   2     3.10738450E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.98586015E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.98476017E+00   # C10 mu+mu-
 03051313 4237   00   2    -8.83778626E-05   # C10' mu+mu-
 03051212 4137   00   0     1.50580593E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50554064E+00   # C11 nu_1 nu_1
 03051212 4237   00   2     2.13400990E-05   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50580593E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50554064E+00   # C11 nu_2 nu_2
 03051414 4237   00   2     2.13401021E-05   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50580593E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50554007E+00   # C11 nu_3 nu_3
 03051616 4237   00   2     2.13403509E-05   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.86059898E-07   # C7
     0305 4422   00   2     2.00799052E-06   # C7
     0305 4322   00   2     1.89229023E-06   # C7'
     0305 6421   00   0     3.30685732E-07   # C8
     0305 6421   00   2    -8.21292605E-06   # C8
     0305 6321   00   2     1.07477917E-06   # C8'
 03051111 4133   00   2    -2.43218990E-07   # C9 e+e-
 03051111 4233   00   2     9.18409074E-08   # C9' e+e-
 03051111 4137   00   2     1.23656855E-05   # C10 e+e-
 03051111 4237   00   2    -2.64541946E-06   # C10' e+e-
 03051313 4133   00   2    -2.43217495E-07   # C9 mu+mu-
 03051313 4233   00   2     9.18408532E-08   # C9' mu+mu-
 03051313 4137   00   2     1.23656871E-05   # C10 mu+mu-
 03051313 4237   00   2    -2.64541964E-06   # C10' mu+mu-
 03051212 4137   00   2    -2.97912769E-06   # C11 nu_1 nu_1
 03051212 4237   00   2     6.38796269E-07   # C11' nu_1 nu_1
 03051414 4137   00   2    -2.97912771E-06   # C11 nu_2 nu_2
 03051414 4237   00   2     6.38796269E-07   # C11' nu_2 nu_2
 03051616 4137   00   2    -2.98950319E-06   # C11 nu_3 nu_3
 03051616 4237   00   2     6.38785053E-07   # C11' nu_3 nu_3
