# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 29.03.2023,  20:08
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.2.3      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
# Either the general MSSM or a model has been used
# which has not yet been implemented in the LesHouches standard
Block MINPAR  # Input parameters
    3    1.68332648E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    7.50000000E+02  # M_1
    2    1.94410000E+03  # M_2
    3    5.83240000E+03  # M_3
   23    2.65310000E+03  # mu
   25    1.63000000E+01  # tan(beta)
   26    1.91790000E+03  # m_A, pole mass
Block SMINPUTS  # SM parameters
         1     1.27908953E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637000E-05  # G_mu [GeV^-2]
         3     1.17200000E-01  # alpha_s(MZ)^MSbar
         4     9.11884000E+01  # m_Z(pole)
         5     4.23000000E+00  # m_b(m_b), MSbar
         6     1.73300000E+02  # m_t(pole)
         7     1.77700000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998910E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658000E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     3.00000000E-03  # m_u(2 GeV), MSbar
        23     1.05000000E-01  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  1.00000000E+03  # (SUSY scale)
   1    3.61013505E-01  # g'(Q)^DRbar
   2    6.32146008E-01  # g(Q)^DRbar
   3    1.01706739E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.57001352E-06   # Y_u(Q)^DRbar
  2  2     3.62797109E-03   # Y_c(Q)^DRbar
  3  3     8.40808489E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.11215350E-04   # Y_d(Q)^DRbar
  2  2     4.43551883E-03   # Y_s(Q)^DRbar
  3  3     2.14245980E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     4.92855106E-05   # Y_e(Q)^DRbar
  2  2     1.01906262E-02   # Y_mu(Q)^DRbar
  3  3     1.71301890E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.16685930E-05   # A_u(Q)^DRbar
  2  2     2.75636155E-08   # A_c(Q)^DRbar
  3  3    -4.38161609E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     4.73450439E-07   # A_d(Q)^DRbar
  2  2     2.25452769E-08   # A_s(Q)^DRbar
  3  3     4.66753215E-10   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.02899389E-06   # A_e(Q)^DRbar
  2  2     9.81293963E-09   # A_mu(Q)^DRbar
  3  3     5.83764720E-10   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    7.50000000E+02  # M_1
   2    1.94410000E+03  # M_2
   3    5.83240000E+03  # M_3
  21   -3.36894889E+06  # M^2_(H,d)
  22   -7.34468212E+06  # M^2_(H,u)
  31    2.24830000E+03  # M_(L,11)
  32    2.24830000E+03  # M_(L,22)
  33    2.24830000E+03  # M_(L,33)
  34    2.24830000E+03  # M_(E,11)
  35    2.24830000E+03  # M_(E,22)
  36    2.24830000E+03  # M_(E,33)
  41    3.05430000E+03  # M_(Q,11)
  42    3.05430000E+03  # M_(Q,22)
  43    2.14370000E+03  # M_(Q,33)
  44    3.05430000E+03  # M_(U,11)
  45    3.05430000E+03  # M_(U,22)
  46    1.97900000E+03  # M_(U,33)
  47    3.05430000E+03  # M_(D,11)
  48    3.05430000E+03  # M_(D,22)
  49    2.14370000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02835389E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.28756499E+02  # h0
        35     1.91738034E+03  # H0
        36     1.91790000E+03  # A0
        37     1.91954571E+03  # H+
   1000001     2.57934975E+03  # ~d_L
   2000001     2.58195855E+03  # ~d_R
   1000002     2.57819652E+03  # ~u_L
   2000002     2.58130395E+03  # ~u_R
   1000003     2.57929326E+03  # ~s_L
   2000003     2.58201612E+03  # ~s_R
   1000004     2.57819660E+03  # ~c_L
   2000004     2.58130516E+03  # ~c_R
   1000005     1.26052678E+03  # ~b_1
   2000005     1.33426947E+03  # ~b_2
   1000006     8.26136944E+02  # ~t_1
   2000006     1.43508660E+03  # ~t_2
   1000011     2.25229191E+03  # ~e_L-
   2000011     2.25039039E+03  # ~e_R-
   1000012     2.25059254E+03  # ~nu_eL
   1000013     2.25275022E+03  # ~mu_L-
   2000013     2.24992893E+03  # ~mu_R-
   1000014     2.25059174E+03  # ~nu_muL
   1000015     2.23339132E+03  # ~tau_1-
   2000015     2.26845546E+03  # ~tau_2-
   1000016     2.25036540E+03  # ~nu_tauL
   1000021     5.28893150E+03  # ~g
   1000022     7.57994356E+02  # ~chi_10
   1000023     1.98669435E+03  # ~chi_20
   1000025    -2.64432563E+03  # ~chi_30
   1000035     2.64935204E+03  # ~chi_40
   1000024     1.98676629E+03  # ~chi_1+
   1000037     2.64942510E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -6.18467103E-02   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    2.65310000E+03  # mu
   2    1.63000000E+01  # tan[beta](Q)
   3    2.45201214E+02  # v(Q)
   4    2.07357825E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     5.15580313E-01   # Re[R_st(1,1)]
   1  2     8.56841258E-01   # Re[R_st(1,2)]
   2  1    -8.56841258E-01   # Re[R_st(2,1)]
   2  2     5.15580313E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     7.10216976E-01   # Re[R_sb(1,1)]
   1  2     7.03982846E-01   # Re[R_sb(1,2)]
   2  1    -7.03982846E-01   # Re[R_sb(2,1)]
   2  2     7.10216976E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.86921370E-01   # Re[R_sta(1,1)]
   1  2     7.26731746E-01   # Re[R_sta(1,2)]
   2  1    -7.26731746E-01   # Re[R_sta(2,1)]
   2  2     6.86921370E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99802313E-01   # Re[N(1,1)]
   1  2    -5.25927401E-04   # Re[N(1,2)]
   1  3     1.87557084E-02   # Re[N(1,3)]
   1  4    -6.57882029E-03   # Re[N(1,4)]
   2  1     2.14764453E-03   # Re[N(2,1)]
   2  2     9.96292067E-01   # Re[N(2,2)]
   2  3    -6.81354764E-02   # Re[N(2,3)]
   2  4     5.24886733E-02   # Re[N(2,4)]
   3  1    -8.60124634E-03   # Re[N(3,1)]
   3  2     1.11028683E-02   # Re[N(3,2)]
   3  3     7.06874183E-01   # Re[N(3,3)]
   3  4     7.07199854E-01   # Re[N(3,4)]
   4  1    -1.77971980E-02   # Re[N(4,1)]
   4  2     8.53145126E-02   # Re[N(4,2)]
   4  3     7.03800163E-01   # Re[N(4,3)]
   4  4    -7.05031932E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.95048718E-01   # Re[U(1,1)]
   1  2     9.93883760E-02   # Re[U(1,2)]
   2  1     9.93883760E-02   # Re[U(2,1)]
   2  2     9.95048718E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.97010093E-01   # Re[V(1,1)]
   1  2     7.72714356E-02   # Re[V(1,2)]
   2  1     7.72714356E-02   # Re[V(2,1)]
   2  2     9.97010093E-01   # Re[V(2,2)]
DECAY   2000011     9.16842855E+00   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99998710E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000011     3.60249925E+00   # ~e^-_L
#    BR                NDA      ID1      ID2
     6.35910445E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     1.21805829E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     2.42283726E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     8.25381114E+00   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.74207936E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     8.62747769E-03    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     1.71645863E-02    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
DECAY   1000013     4.51706894E+00   # ~mu^-_L
#    BR                NDA      ID1      ID2
     7.56763296E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     8.13773746E-02    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     1.61859329E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     6.45964817E+00   # ~tau^-_1
#    BR                NDA      ID1      ID2
     9.11147450E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     2.96938060E-02    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     5.91587436E-02    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
DECAY   2000015     6.30072609E+00   # ~tau^-_2
#    BR                NDA      ID1      ID2
     8.84041712E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     3.88495228E-02    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     7.71070387E-02    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     3.59472988E+00   # ~nu_e
#    BR                NDA      ID1      ID2
     6.38894365E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     1.20096405E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     2.41009230E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000014     3.59472276E+00   # ~nu_mu
#    BR                NDA      ID1      ID2
     6.38895286E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     1.20095996E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     2.41008718E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000016     3.59271357E+00   # ~nu_tau
#    BR                NDA      ID1      ID2
     6.39155390E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     1.19980505E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     2.40863817E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000001     5.35211412E+00   # ~d_L
#    BR                NDA      ID1      ID2
     5.82547496E-02    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     3.14310295E-01    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     6.27434955E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
DECAY   2000001     1.24202725E+00   # ~d_R
#    BR                NDA      ID1      ID2
     9.99793935E-01    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     1.36663795E-04    2    -1000024         2   # BR(~d_R -> chi^-_1 u)
DECAY   1000003     5.26665285E+00   # ~s_L
#    BR                NDA      ID1      ID2
     6.29310633E-02    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     3.12749173E-01    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     6.24319764E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
DECAY   2000003     1.32747640E+00   # ~s_R
#    BR                NDA      ID1      ID2
     9.20635552E-01    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     2.64941690E-02    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     5.28702786E-02    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
DECAY   1000005     2.93137653E+01   # ~b_1
#    BR                NDA      ID1      ID2
     6.61124002E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     9.93388760E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     4.36849995E+01   # ~b_2
#    BR                NDA      ID1      ID2
     4.81357459E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     9.95186425E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   1000002     5.34927204E+00   # ~u_L
#    BR                NDA      ID1      ID2
     5.76006820E-02    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     3.14047538E-01    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     6.28351780E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
DECAY   2000002     4.96554617E+00   # ~u_R
#    BR                NDA      ID1      ID2
     9.99999081E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
DECAY   1000004     5.34911756E+00   # ~c_L
#    BR                NDA      ID1      ID2
     5.77574664E-02    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     3.13994411E-01    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     6.28248122E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
DECAY   2000004     4.96570289E+00   # ~c_R
#    BR                NDA      ID1      ID2
     9.99800145E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     1.32738215E-04    2     1000024         3   # BR(~c_R -> chi^+_1 s)
DECAY   1000006     1.74033490E-09   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     7.36578946E+01   # ~t_2
#    BR                NDA      ID1      ID2
     5.74419077E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     7.64221830E-02    2     1000005        24   # BR(~t_2 -> ~b_1 W^+)
     4.95814390E-03    2     2000005        24   # BR(~t_2 -> ~b_2 W^+)
     7.15314346E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     1.97561136E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     1.75566799E+01   # chi^+_1
#    BR                NDA      ID1      ID2
     3.35043056E-01    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     1.98555072E-01    2     2000006        -5   # BR(chi^+_1 -> ~t_2 b_bar)
     2.69447490E-01    2    -1000005         6   # BR(chi^+_1 -> ~b^*_1 t)
     1.95983039E-01    2    -2000005         6   # BR(chi^+_1 -> ~b^*_2 t)
     7.14058678E-04    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     1.62996534E-04    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
DECAY   1000037     8.89763293E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     1.61618396E-04    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     5.58901728E-04    2    -2000015        16   # BR(chi^+_2 -> ~tau^+_2 nu_tau)
     7.20738955E-04    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     3.57032014E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     1.14162247E-01    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     2.04982364E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     1.70742937E-01    2    -2000005         6   # BR(chi^+_2 -> ~b^*_2 t)
     1.93238317E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     4.45349343E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     4.24733815E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     4.29681710E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     4.37556413E-04    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     1.00293415E-03    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     5.63718571E-04    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     1.75010534E+01   # chi^0_2
#    BR                NDA      ID1      ID2
     1.72014573E-01    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     1.72014573E-01    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     1.04600996E-01    2     2000006        -6   # BR(chi^0_2 -> ~t_2 t_bar)
     1.04600996E-01    2    -2000006         6   # BR(chi^0_2 -> ~t^*_2 t)
     1.29018733E-01    2     1000005        -5   # BR(chi^0_2 -> ~b_1 b_bar)
     1.29018733E-01    2    -1000005         5   # BR(chi^0_2 -> ~b^*_1 b)
     9.38923922E-02    2     2000005        -5   # BR(chi^0_2 -> ~b_2 b_bar)
     9.38923922E-02    2    -2000005         5   # BR(chi^0_2 -> ~b^*_2 b)
     5.84526923E-04    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.89043932E-04    3     1000022         6        -6   # BR(chi^0_2 -> chi^0_1 t t_bar)
DECAY   1000025     8.88694828E+01   # chi^0_3
#    BR                NDA      ID1      ID2
     3.34320172E-04    2     1000015       -15   # BR(chi^0_3 -> ~tau^-_1 tau^+)
     3.34320172E-04    2    -1000015        15   # BR(chi^0_3 -> ~tau^+_1 tau^-)
     3.22411992E-04    2     2000015       -15   # BR(chi^0_3 -> ~tau^-_2 tau^+)
     3.22411992E-04    2    -2000015        15   # BR(chi^0_3 -> ~tau^+_2 tau^-)
     2.25081235E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     2.25081235E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     1.78780034E-01    2     2000006        -6   # BR(chi^0_3 -> ~t_2 t_bar)
     1.78780034E-01    2    -2000006         6   # BR(chi^0_3 -> ~t^*_2 t)
     1.15479474E-02    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     1.15479474E-02    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     1.19096251E-02    2     2000005        -5   # BR(chi^0_3 -> ~b_2 b_bar)
     1.19096251E-02    2    -2000005         5   # BR(chi^0_3 -> ~b^*_2 b)
     4.11167273E-02    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     4.11167273E-02    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     1.63175945E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     3.85149201E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     3.92762556E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     6.53764185E-04    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     9.26648459E-04    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     2.94390579E-04    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     5.58292657E-04    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     5.58292657E-04    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     8.89926664E+01   # chi^0_4
#    BR                NDA      ID1      ID2
     2.48433735E-04    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     2.48433735E-04    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     4.26153264E-04    2     2000015       -15   # BR(chi^0_4 -> ~tau^-_2 tau^+)
     4.26153264E-04    2    -2000015        15   # BR(chi^0_4 -> ~tau^+_2 tau^-)
     2.68926896E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     2.68926896E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     1.37394048E-01    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     1.37394048E-01    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     9.12751321E-03    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     9.12751321E-03    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     1.47626971E-02    2     2000005        -5   # BR(chi^0_4 -> ~b_2 b_bar)
     1.47626971E-02    2    -2000005         5   # BR(chi^0_4 -> ~b^*_2 b)
     3.61872075E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     3.61872075E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     3.76708727E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     6.83802996E-04    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     1.62557936E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     4.23279172E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.06082771E-03    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     1.66389078E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     5.71224731E-04    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     5.71224731E-04    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     8.95104421E+02   # ~g
#    BR                NDA      ID1      ID2
     3.53369154E-02    2     1000002        -2   # BR(~g -> ~u_L u_bar)
     3.53369154E-02    2    -1000002         2   # BR(~g -> ~u^*_L u)
     3.52838023E-02    2     2000002        -2   # BR(~g -> ~u_R u_bar)
     3.52838023E-02    2    -2000002         2   # BR(~g -> ~u^*_R u)
     3.53375134E-02    2     1000004        -4   # BR(~g -> ~c_L c_bar)
     3.53375134E-02    2    -1000004         4   # BR(~g -> ~c^*_L c)
     3.52831795E-02    2     2000004        -4   # BR(~g -> ~c_R c_bar)
     3.52831795E-02    2    -2000004         4   # BR(~g -> ~c^*_R c)
     5.44342914E-02    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     5.44342914E-02    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     5.54227470E-02    2     2000006        -6   # BR(~g -> ~t_2 t_bar)
     5.54227470E-02    2    -2000006         6   # BR(~g -> ~t^*_2 t)
     3.53172056E-02    2     1000001        -1   # BR(~g -> ~d_L d_bar)
     3.53172056E-02    2    -1000001         1   # BR(~g -> ~d^*_L d)
     3.52726120E-02    2     2000001        -1   # BR(~g -> ~d_R d_bar)
     3.52726120E-02    2    -2000001         1   # BR(~g -> ~d^*_R d)
     3.53176453E-02    2     1000003        -3   # BR(~g -> ~s_L s_bar)
     3.53176453E-02    2    -1000003         3   # BR(~g -> ~s^*_L s)
     3.52721530E-02    2     2000003        -3   # BR(~g -> ~s_R s_bar)
     3.52721530E-02    2    -2000003         3   # BR(~g -> ~s^*_R s)
     5.39960366E-02    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     5.39960366E-02    2    -1000005         5   # BR(~g -> ~b^*_1 b)
     5.33971273E-02    2     2000005        -5   # BR(~g -> ~b_2 b_bar)
     5.33971273E-02    2    -2000005         5   # BR(~g -> ~b^*_2 b)
#    BR                NDA      ID1      ID2       ID3
     1.64311315E-04    3     1000022         6        -6   # BR(~g -> chi^0_1 t t_bar)
DECAY        25     3.66503589E-03   # h^0
#    BR                NDA      ID1      ID2
     2.87375523E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     8.09222587E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.09478766E-04    2           3        -3   # BR(h^0 -> s s_bar)
     4.76475823E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.47092642E-02    2           4        -4   # BR(h^0 -> c c_bar)
     9.05531978E-02    2          21        21   # BR(h^0 -> g g)
     3.41623306E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.39929551E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.39929551E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.39929551E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     4.89753428E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     4.89753428E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.39929551E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.39929551E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.39929551E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     4.89753428E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     4.89753428E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     1.17483072E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     1.17483072E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     1.00699776E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     6.71331842E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     4.69932290E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     4.69932290E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     4.69932290E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     4.69932290E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     4.69932290E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     3.24577647E+01   # H^0
#    BR                NDA      ID1      ID2
     1.21579206E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     3.43541840E-02    2          15       -15   # BR(H^0 -> tau^- tau^+)
     1.83766918E-01    2           5        -5   # BR(H^0 -> b b_bar)
     1.05786157E-02    2           6        -6   # BR(H^0 -> t t_bar)
     7.70040944E-01    2     1000006  -1000006   # BR(H^0 -> ~t_1 ~t^*_1)
     2.14931492E-04    2          23        23   # BR(H^0 -> Z Z)
     4.08060928E-04    2          24       -24   # BR(H^0 -> W^+ W^-)
     3.97256059E-04    2          21        21   # BR(H^0 -> g g)
DECAY        36     7.43474987E+00   # A^0
#    BR                NDA      ID1      ID2
     5.30957748E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.50031270E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     3.43978672E-04    2           3        -3   # BR(A^0 -> s s_bar)
     8.02544710E-01    2           5        -5   # BR(A^0 -> b b_bar)
     4.62793544E-02    2           6        -6   # BR(A^0 -> t t_bar)
     1.50046329E-04    2     1000022   1000022   # BR(A^0 -> chi^0_1 chi^0_1)
DECAY        37     6.57491703E+00   # H^+
#    BR                NDA      ID1      ID2
     6.00908775E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.69797140E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     3.42380266E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.29151604E-01    2          -5         6   # BR(H^+ -> b_bar t)
     1.07179042E-04    2          25        24   # BR(H^+ -> h^0 W^+)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.01879344E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    2.65671207E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    2.65690000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99929266E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    3.83451933E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    3.76378486E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.01879344E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    2.65671207E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    2.65690000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999671E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    3.29260971E-07        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999671E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    3.29260971E-07        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15101023E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    2.62384165E-01        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.04067246E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    3.29260971E-07        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999671E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.01898012E-04   # BR(b -> s gamma)
    2    1.59205485E-06   # BR(b -> s mu+ mu-)
    3    3.54402400E-05   # BR(b -> s nu nu)
    4    5.19095457E-11   # BR(Bd -> mu+ mu-)
    5    3.39130149E-09   # BR(Bs -> mu+ mu-)
    6    1.10286952E-04   # BR(B_u -> tau nu)
    7    9.97263943E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73951757E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.03304044E+01   # |Delta(M_Bs)| [ps^-1] 
   20    1.09204091E-14   # Delta(g-2)_electron/2
   21    4.08502708E-11   # Delta(g-2)_muon/2
   22   -4.92477056E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    4.60869149E-03   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87801604E-01   # C7
     0305 4422   00   1    -3.02372924E-03   # C7
     0305 4322   00   1    -6.56247558E-05   # C7'
     0305 6421   00   0    -9.49035011E-02   # C8
     0305 6421   00   1    -7.69014470E-03   # C8
     0305 6321   00   1    -1.62221853E-04   # C8'
 03051111 4133   00   0    -4.84461120E-01   # C9 e+e-
 03051111 4133   00   1     1.60851571E-04   # C9 e+e-
 03051111 4233   00   1    -3.51190140E-06   # C9' e+e-
 03051111 4137   00   0    -3.92552721E+00   # C10 e+e-
 03051111 4137   00   1     1.09100349E-03   # C10 e+e-
 03051111 4237   00   1     1.07080487E-04   # C10' e+e-
 03051313 4133   00   0    -4.84461120E-01   # C9 mu+mu-
 03051313 4133   00   1     1.60841748E-04   # C9 mu+mu-
 03051313 4233   00   1    -3.51194728E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.92552721E+00   # C10 mu+mu-
 03051313 4137   00   1     1.09099366E-03   # C10 mu+mu-
 03051313 4237   00   1     1.07080533E-04   # C10' mu+mu-
 03051212 4237   00   0     1.49077528E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -2.64587855E-04   # C11 nu_1 nu_1
 03051212 4137   00   1    -2.58953913E-05   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49077528E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -2.64587855E-04   # C11 nu_2 nu_2
 03051414 4137   00   1    -2.58953805E-05   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49077528E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -2.64587846E-04   # C11 nu_3 nu_3
 03051616 4137   00   1    -2.58923353E-05   # C11' nu_3 nu_3
