# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 09.06.2023,  11:23
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
    3    1.38621799E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    7.31000000E+02  # M_1
    2    7.66000000E+02  # M_2
    3    1.90630000E+03  # M_3
   23    1.28630000E+03  # mu
   25    1.34000000E+01  # tan(beta)
   26    1.59290000E+03  # m_A, pole mass
   31    3.58960000E+03  # M^2_L11
   32    3.58960000E+03  # M^2_L22
   33    3.58960000E+03  # M^2_L33
   34    3.58960000E+03  # M^2_E11
   35    3.58960000E+03  # M^2_E22
   36    3.58960000E+03  # M^2_E33
   41    3.25260000E+03  # M^2_Q11
   42    3.25260000E+03  # M^2_Q22
   43    1.63430000E+03  # M^2_Q33
   44    3.25260000E+03  # M^2_U11
   45    3.25260000E+03  # M^2_U22
   46    1.05440000E+03  # M^2_U33
   47    3.25260000E+03  # M^2_D11
   48    3.25260000E+03  # M^2_D22
   49    1.63430000E+03  # M^2_D33
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
   1    3.61188547E-01  # g'(Q)^DRbar
   2    6.35312108E-01  # g(Q)^DRbar
   3    1.03727668E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.65073699E-06   # Y_u(Q)^DRbar
  2  2     3.66214424E-03   # Y_c(Q)^DRbar
  3  3     8.51065329E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.86547821E-04   # Y_d(Q)^DRbar
  2  2     3.91750111E-03   # Y_s(Q)^DRbar
  3  3     1.87186141E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     3.98960519E-05   # Y_e(Q)^DRbar
  2  2     8.24920001E-03   # Y_mu(Q)^DRbar
  3  3     1.38694002E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3    -3.28094672E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     0.00000000E+00   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     0.00000000E+00   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    7.31000000E+02  # M_1
   2    7.66000000E+02  # M_2
   3    1.90630000E+03  # M_3
  21    8.37497460E+05  # M^2_(H,d)
  22   -1.73717329E+06  # M^2_(H,u)
  31    3.58960000E+03  # M_(L,11)
  32    3.58960000E+03  # M_(L,22)
  33    3.58960000E+03  # M_(L,33)
  34    3.58960000E+03  # M_(E,11)
  35    3.58960000E+03  # M_(E,22)
  36    3.58960000E+03  # M_(E,33)
  41    3.25260000E+03  # M_(Q,11)
  42    3.25260000E+03  # M_(Q,22)
  43    1.63430000E+03  # M_(Q,33)
  44    3.25260000E+03  # M_(U,11)
  45    3.25260000E+03  # M_(U,22)
  46    1.05440000E+03  # M_(U,33)
  47    3.25260000E+03  # M_(D,11)
  48    3.25260000E+03  # M_(D,22)
  49    1.63430000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02829600E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.26264062E+02  # h0
        35     1.59266396E+03  # H0
        36     1.59290000E+03  # A0
        37     1.59510083E+03  # H+
   1000001     3.27161821E+03  # ~d_L
   2000001     3.26583698E+03  # ~d_R
   1000002     3.27086933E+03  # ~u_L
   2000002     3.26621254E+03  # ~u_R
   1000003     3.27162151E+03  # ~s_L
   2000003     3.26583472E+03  # ~s_R
   1000004     3.27086994E+03  # ~c_L
   2000004     3.26621309E+03  # ~c_R
   1000005     1.63100406E+03  # ~b_1
   2000005     1.66821717E+03  # ~b_2
   1000006     1.00903220E+03  # ~t_1
   2000006     1.68511087E+03  # ~t_2
   1000011     3.59585985E+03  # ~e_L-
   2000011     3.59263235E+03  # ~e_R-
   1000012     3.59463332E+03  # ~nu_eL
   1000013     3.59588128E+03  # ~mu_L-
   2000013     3.59261525E+03  # ~mu_R-
   1000014     3.59463477E+03  # ~nu_muL
   1000015     3.59035081E+03  # ~tau_1-
   2000015     3.59936745E+03  # ~tau_2-
   1000016     3.59504530E+03  # ~nu_tauL
   1000021     2.04994572E+03  # ~g
   1000022     7.38116743E+02  # ~chi_10
   1000023     8.02389529E+02  # ~chi_20
   1000025    -1.28843447E+03  # ~chi_30
   1000035     1.29450637E+03  # ~chi_40
   1000024     8.02336006E+02  # ~chi_1+
   1000037     1.29507533E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -7.53431908E-02   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    1.28630000E+03  # mu
   2    1.34000000E+01  # tan[beta](Q)
   3    2.43737301E+02  # v(Q)
   4    2.51240825E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     2.87029205E-01   # Re[R_st(1,1)]
   1  2     9.57921832E-01   # Re[R_st(1,2)]
   2  1    -9.57921832E-01   # Re[R_st(2,1)]
   2  2     2.87029205E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.42847237E-01   # Re[R_sb(1,1)]
   1  2     3.33225282E-01   # Re[R_sb(1,2)]
   2  1    -3.33225282E-01   # Re[R_sb(2,1)]
   2  2     9.42847237E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     5.88915387E-01   # Re[R_sta(1,1)]
   1  2     8.08194696E-01   # Re[R_sta(1,2)]
   2  1    -8.08194696E-01   # Re[R_sta(2,1)]
   2  2     5.88915387E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.96329795E-01   # Re[N(1,1)]
   1  2     4.90757344E-02   # Re[N(1,2)]
   1  3    -5.94074776E-02   # Re[N(1,3)]
   1  4     3.72728270E-02   # Re[N(1,4)]
   2  1    -5.71952662E-02   # Re[N(2,1)]
   2  2    -9.91252111E-01   # Re[N(2,2)]
   2  3     9.90386269E-02   # Re[N(2,3)]
   2  4    -6.58733917E-02   # Re[N(2,4)]
   3  1    -1.44357834E-02   # Re[N(3,1)]
   3  2     2.43670896E-02   # Re[N(3,2)]
   3  3     7.06133489E-01   # Re[N(3,3)]
   3  4     7.07512083E-01   # Re[N(3,4)]
   4  1    -6.20261994E-02   # Re[N(4,1)]
   4  2     1.20071100E-01   # Re[N(4,2)]
   4  3     6.98596877E-01   # Re[N(4,3)]
   4  4    -7.02636524E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.89437892E-01   # Re[U(1,1)]
   1  2     1.44957432E-01   # Re[U(1,2)]
   2  1     1.44957432E-01   # Re[U(2,1)]
   2  2     9.89437892E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.95315051E-01   # Re[V(1,1)]
   1  2     9.66847982E-02   # Re[V(1,2)]
   2  1     9.66847982E-02   # Re[V(2,1)]
   2  2     9.95315051E-01   # Re[V(2,2)]
DECAY   2000011     1.70942731E+01   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.93427528E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
     3.22167621E-03    2     1000023        11   # BR(~e^-_R -> chi^0_2 e^-)
     1.72620626E-04    2     1000025        11   # BR(~e^-_R -> chi^0_3 e^-)
     3.17794350E-03    2     1000035        11   # BR(~e^-_R -> chi^0_4 e^-)
DECAY   1000011     4.32713323E+01   # ~e^-_L
#    BR                NDA      ID1      ID2
     8.19350647E-02    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     3.15731165E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     1.81797682E-03    2     1000035        11   # BR(~e^-_L -> chi^0_4 e^-)
     5.89829890E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
     1.06197094E-02    2    -1000037        12   # BR(~e^-_L -> chi^-_2 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.72658793E+01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.78796912E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     8.28751912E-03    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     2.35419278E-04    2     1000025        13   # BR(~mu^-_R -> chi^0_3 mu^-)
     3.06111441E-03    2     1000035        13   # BR(~mu^-_R -> chi^0_4 mu^-)
     9.61661257E-03    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
DECAY   1000013     4.31108590E+01   # ~mu^-_L
#    BR                NDA      ID1      ID2
     8.41463249E-02    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     3.14866749E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     1.26109107E-04    2     1000025        13   # BR(~mu^-_L -> chi^0_3 mu^-)
     1.94246561E-03    2     1000035        13   # BR(~mu^-_L -> chi^0_4 mu^-)
     5.88176255E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     1.07420962E-02    2    -1000037        14   # BR(~mu^-_L -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     2.80731969E+01   # ~tau^-_1
#    BR                NDA      ID1      ID2
     4.43818406E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.83122715E-01    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     1.59986808E-02    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     7.91403813E-03    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     3.42962433E-01    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
     6.18372685E-03    2    -1000037        16   # BR(~tau^-_1 -> chi^-_2 nu_tau)
DECAY   2000015     3.54409473E+01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     2.28027715E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     2.42690480E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     1.69016352E-02    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     2.62317291E-02    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     4.49245047E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     3.69033935E-02    2    -1000037        16   # BR(~tau^-_2 -> chi^-_2 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     4.32887452E+01   # ~nu_e
#    BR                NDA      ID1      ID2
     1.15879719E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     2.76666928E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     2.68695631E-04    2     1000025        12   # BR(~nu_e -> chi^0_3 nu_e)
     6.09310966E-03    2     1000035        12   # BR(~nu_e -> chi^0_4 nu_e)
     5.96371603E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
     4.71994449E-03    2     1000037        11   # BR(~nu_e -> chi^+_2 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000014     4.32924662E+01   # ~nu_mu
#    BR                NDA      ID1      ID2
     1.15869814E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     2.76643284E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     2.68672709E-04    2     1000025        14   # BR(~nu_mu -> chi^0_3 nu_mu)
     6.09258989E-03    2     1000035        14   # BR(~nu_mu -> chi^0_4 nu_mu)
     5.96322737E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     4.80290204E-03    2     1000037        13   # BR(~nu_mu -> chi^+_2 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000016     4.43408103E+01   # ~nu_tau
#    BR                NDA      ID1      ID2
     1.13145514E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     2.70139959E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     2.62368137E-04    2     1000025        16   # BR(~nu_tau -> chi^0_3 nu_tau)
     5.94962783E-03    2     1000035        16   # BR(~nu_tau -> chi^0_4 nu_tau)
     5.82882068E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     2.76204629E-02    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000001     7.01532092E+01   # ~d_R
#    BR                NDA      ID1      ID2
     2.39991951E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.75845289E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     1.03996826E+02   # ~d_L
#    BR                NDA      ID1      ID2
     6.43847935E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     1.07238464E-01    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.56144748E-03    2     1000035         1   # BR(~d_L -> chi^0_4 d)
     2.18447346E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     3.77480931E-03    2    -1000037         2   # BR(~d_L -> chi^-_2 u)
     6.62473213E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     7.01711435E+01   # ~s_R
#    BR                NDA      ID1      ID2
     2.39884061E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     1.60863616E-04    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     1.70438241E-04    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
     9.75596308E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     1.03981730E+02   # ~s_L
#    BR                NDA      ID1      ID2
     6.44261642E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     1.07198096E-01    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.56801828E-03    2     1000035         3   # BR(~s_L -> chi^0_4 s)
     2.18364446E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     3.78654691E-03    2    -1000037         4   # BR(~s_L -> chi^-_2 c)
     6.62569850E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     5.78741671E+01   # ~b_1
#    BR                NDA      ID1      ID2
     5.44623536E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     5.72621686E-02    2     1000023         5   # BR(~b_1 -> chi^0_2 b)
     1.30553682E-03    2     1000025         5   # BR(~b_1 -> chi^0_3 b)
     1.11542532E-03    2     1000035         5   # BR(~b_1 -> chi^0_4 b)
     1.09889862E-01    2    -1000024         6   # BR(~b_1 -> chi^-_1 t)
     4.17196255E-02    2    -1000037         6   # BR(~b_1 -> chi^-_2 t)
     7.83261147E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     9.07701304E+00   # ~b_2
#    BR                NDA      ID1      ID2
     5.92714386E-02    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     3.65888387E-02    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     1.10663413E-02    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     1.35842835E-02    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     6.91900214E-02    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     8.44906255E-02    2    -1000037         6   # BR(~b_2 -> chi^-_2 t)
     7.25808451E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   2000002     7.52659325E+01   # ~u_R
#    BR                NDA      ID1      ID2
     8.94883799E-02    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     2.89163464E-04    2     1000023         2   # BR(~u_R -> chi^0_2 u)
     2.73662423E-04    2     1000035         2   # BR(~u_R -> chi^0_4 u)
     9.09933919E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     1.03987778E+02   # ~u_L
#    BR                NDA      ID1      ID2
     2.22112510E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     1.12010839E-01    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.05385237E-03    2     1000035         2   # BR(~u_L -> chi^0_4 u)
     2.21006103E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     1.67879418E-03    2     1000037         1   # BR(~u_L -> chi^+_2 d)
     6.61987110E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     7.52676215E+01   # ~c_R
#    BR                NDA      ID1      ID2
     8.94860009E-02    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     2.90523599E-04    2     1000023         4   # BR(~c_R -> chi^0_2 c)
     2.77208075E-04    2     1000035         4   # BR(~c_R -> chi^0_4 c)
     9.09916890E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     1.03988647E+02   # ~c_L
#    BR                NDA      ID1      ID2
     2.22137979E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     1.12008942E-01    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.05721840E-03    2     1000035         4   # BR(~c_L -> chi^0_4 c)
     2.21002545E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     1.68593498E-03    2     1000037         3   # BR(~c_L -> chi^+_2 s)
     6.61978881E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     5.62756410E-01   # ~t_1
#    BR                NDA      ID1      ID2
     5.90969296E-01    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     8.22836942E-02    2     1000023         6   # BR(~t_1 -> chi^0_2 t)
     3.26736597E-01    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000006     6.57649643E+01   # ~t_2
#    BR                NDA      ID1      ID2
     2.65299515E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     4.98549403E-02    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     2.19151070E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     3.60305136E-02    2     1000035         6   # BR(~t_2 -> chi^0_4 t)
     1.03467492E-01    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     1.05632782E-02    2     1000037         5   # BR(~t_2 -> chi^+_2 b)
     4.50375889E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     3.25139784E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     2.16136542E-05   # chi^+_1
#    BR                NDA      ID1      ID2       ID3
     3.33680700E-01    3     1000022        -1         2   # BR(chi^+_1 -> chi^0_1 d_bar u)
     3.33272986E-01    3     1000022        -3         4   # BR(chi^+_1 -> chi^0_1 s_bar c)
     1.11148640E-01    3     1000022       -11        12   # BR(chi^+_1 -> chi^0_1 e^+ nu_e)
     1.11148117E-01    3     1000022       -13        14   # BR(chi^+_1 -> chi^0_1 mu^+ nu_mu)
     1.10749556E-01    3     1000022       -15        16   # BR(chi^+_1 -> chi^0_1 tau^+ nu_tau)
DECAY   1000037     1.15144092E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     3.28744225E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     5.80801300E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     2.11831040E-01    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.00893614E-01    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     1.97818948E-01    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     6.26620989E-04    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     5.24068856E-04    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     1.47605957E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     6.82420148E-08   # chi^0_2
#    BR                NDA      ID1      ID2
     3.59492695E-03    2     1000022        22   # BR(chi^0_2 -> chi^0_1 photon)
#    BR                NDA      ID1      ID2       ID3
     1.15403855E-01    3     1000022         2        -2   # BR(chi^0_2 -> chi^0_1 u u_bar)
     1.15586553E-01    3     1000022         4        -4   # BR(chi^0_2 -> chi^0_1 c c_bar)
     1.52711155E-01    3     1000022         1        -1   # BR(chi^0_2 -> chi^0_1 d d_bar)
     1.52713038E-01    3     1000022         3        -3   # BR(chi^0_2 -> chi^0_1 s s_bar)
     2.00005364E-01    3     1000022         5        -5   # BR(chi^0_2 -> chi^0_1 b b_bar)
     2.86454426E-02    3     1000022        11       -11   # BR(chi^0_2 -> chi^0_1 e^- e^+)
     2.86519227E-02    3     1000022        13       -13   # BR(chi^0_2 -> chi^0_1 mu^- mu^+)
     3.04520351E-02    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
     1.72235708E-01    3     1000022        12       -12   # BR(chi^0_2 -> chi^0_1 nu_e nu_bar_e)
DECAY   1000025     9.77412972E+00   # chi^0_3
#    BR                NDA      ID1      ID2
     1.11088428E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     1.11088428E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     2.29393869E-01    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     2.29393869E-01    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     9.51756203E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     2.06995332E-01    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     4.84472682E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     8.16735954E-03    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.30990913E-04    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     4.05935893E-04    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     1.55491438E-03    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     1.55491438E-03    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     1.19545975E+01   # chi^0_4
#    BR                NDA      ID1      ID2
     1.83761565E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     1.83761565E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     1.86980683E-01    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     1.86980683E-01    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     4.32969006E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     7.43521610E-03    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     7.31926689E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     1.70576333E-01    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.31022766E-04    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     1.32709458E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     1.35672754E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     1.35672754E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     3.82257895E+01   # ~g
#    BR                NDA      ID1      ID2
     2.87657173E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     2.87657173E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     6.94319063E-02    2     2000006        -6   # BR(~g -> ~t_2 t_bar)
     6.94319063E-02    2    -2000006         6   # BR(~g -> ~t^*_2 t)
     7.67424202E-02    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     7.67424202E-02    2    -1000005         5   # BR(~g -> ~b^*_1 b)
     6.59738947E-02    2     2000005        -5   # BR(~g -> ~b_2 b_bar)
     6.59738947E-02    2    -2000005         5   # BR(~g -> ~b^*_2 b)
#    BR                NDA      ID1      ID2       ID3
DECAY        25     3.58373389E-03   # h^0
#    BR                NDA      ID1      ID2
     2.80619783E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     7.91090540E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.45910962E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.47431136E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.56986752E-02    2           4        -4   # BR(h^0 -> c c_bar)
     8.66154780E-02    2          21        21   # BR(h^0 -> g g)
     3.23393907E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.11106218E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.11106218E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.11106218E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     3.88871761E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     3.88871761E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.11106218E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.11106218E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.11106218E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     3.88871761E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     3.88871761E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     8.81019594E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     8.81019594E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     7.55159652E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     5.03439768E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     3.52407838E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     3.52407838E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     3.52407838E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     3.52407838E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     3.52407838E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     4.84203332E+00   # H^0
#    BR                NDA      ID1      ID2
     4.42774427E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     1.25161566E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     3.41929677E-04    2           3        -3   # BR(H^0 -> s s_bar)
     7.80669677E-01    2           5        -5   # BR(H^0 -> b b_bar)
     8.83972760E-02    2           6        -6   # BR(H^0 -> t t_bar)
     1.68478510E-04    2     1000022   1000022   # BR(H^0 -> chi^0_1 chi^0_1)
     2.69888900E-04    2     1000022   1000023   # BR(H^0 -> chi^0_1 chi^0_2)
     1.44842495E-03    2          23        23   # BR(H^0 -> Z Z)
     2.76598502E-03    2          24       -24   # BR(H^0 -> W^+ W^-)
     2.76415527E-04    2          25        25   # BR(H^0 -> h^0 h^0)
DECAY        36     4.85559735E+00   # A^0
#    BR                NDA      ID1      ID2
     4.41659632E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.24847060E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     3.41069325E-04    2           3        -3   # BR(A^0 -> s s_bar)
     7.78705836E-01    2           5        -5   # BR(A^0 -> b b_bar)
     8.89160138E-02    2           6        -6   # BR(A^0 -> t t_bar)
     1.44188100E-03    2     1000022   1000022   # BR(A^0 -> chi^0_1 chi^0_1)
     5.07807545E-03    2     1000022   1000023   # BR(A^0 -> chi^0_1 chi^0_2)
     1.88569415E-04    2          25        23   # BR(A^0 -> h^0 Z)
DECAY        37     4.24644119E+00   # H^+
#    BR                NDA      ID1      ID2
     5.05713896E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.42953734E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     3.43818046E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.50353135E-01    2          -5         6   # BR(H^+ -> b_bar t)
     5.63121807E-03    2     1000024   1000022   # BR(H^+ -> chi^+_1 chi^0_1)
     2.11593772E-04    2          25        24   # BR(H^+ -> h^0 W^+)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.02302827E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    1.79536972E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    1.79560000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99871752E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    5.69741742E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    5.56916908E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.02302827E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    1.79536972E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    1.79560000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999270E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    7.30000141E-07        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999270E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    7.30000141E-07        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.13090718E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    5.28396345E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.09881377E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    7.30000141E-07        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999270E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    2.84219019E-04   # BR(b -> s gamma)
    2    1.59158658E-06   # BR(b -> s mu+ mu-)
    3    3.54492859E-05   # BR(b -> s nu nu)
    4    5.19220889E-11   # BR(Bd -> mu+ mu-)
    5    3.39196504E-09   # BR(Bs -> mu+ mu-)
    6    1.10264264E-04   # BR(B_u -> tau nu)
    7    9.97058784E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73462192E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.02923711E+01   # |Delta(M_Bs)| [ps^-1] 
   20    1.54035587E-14   # Delta(g-2)_electron/2
   21    3.33599140E-11   # Delta(g-2)_muon/2
   22   -2.08070869E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    1.15383112E-02   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87802455E-01   # C7
     0305 4422   00   1     1.88656213E-02   # C7
     0305 4322   00   1     3.84543419E-04   # C7'
     0305 6421   00   0    -9.49037728E-02   # C8
     0305 6421   00   1     1.82243661E-03   # C8
     0305 6321   00   1     2.98415441E-05   # C8'
 03051111 4133   00   0    -4.57850325E-01   # C9 e+e-
 03051111 4133   00   1     3.39134655E-05   # C9 e+e-
 03051111 4233   00   1    -4.00564728E-06   # C9' e+e-
 03051111 4137   00   0    -3.95218028E+00   # C10 e+e-
 03051111 4137   00   1     3.44320969E-04   # C10 e+e-
 03051111 4237   00   1     1.00382662E-04   # C10' e+e-
 03051313 4133   00   0    -4.57850325E-01   # C9 mu+mu-
 03051313 4133   00   1     3.39095981E-05   # C9 mu+mu-
 03051313 4233   00   1    -4.00568950E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.95218028E+00   # C10 mu+mu-
 03051313 4137   00   1     3.44317024E-04   # C10 mu+mu-
 03051313 4237   00   1     1.00382704E-04   # C10' mu+mu-
 03051212 4237   00   0     1.49078738E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -8.64143300E-05   # C11 nu_1 nu_1
 03051212 4137   00   1    -2.41122171E-05   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49078738E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -8.64143300E-05   # C11 nu_2 nu_2
 03051414 4137   00   1    -2.41122075E-05   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49078738E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -8.64143289E-05   # C11 nu_3 nu_3
 03051616 4137   00   1    -2.41095270E-05   # C11' nu_3 nu_3
