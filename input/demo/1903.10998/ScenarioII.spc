# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 29.03.2023,  12:40
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
    3    5.37752302E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.62919000E+03  # M_1
    2    3.61340000E+03  # M_2
    3    1.72078000E+03  # M_3
   11   -4.43488000E+03  # A_t
   12    2.20173000E+03  # A_b
   13    3.64185000E+03  # A_tau
   23    2.61534000E+03  # mu
   25    5.31300000E+01  # tan(beta)
   26    3.45134000E+03  # m_A, pole mass
Block SMINPUTS  # SM parameters
         1     1.27908970E+02  # alpha_em^-1(MZ)^MSbar
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
Block gauge Q=  2.44796300E+03  # (SUSY scale)
   1    3.63803300E-01  # g'(Q)^DRbar
   2    6.34378008E-01  # g(Q)^DRbar
   3    1.02254260E+00  # g3(Q)^DRbar
Block Yu Q=  2.44796300E+03  # (SUSY scale)
  1  1     6.85587569E-06   # Y_u(Q)^DRbar
  2  2     3.48268081E-03   # Y_c(Q)^DRbar
  3  3     8.30273745E-01   # Y_t(Q)^DRbar
Block Yd Q=  2.44796300E+03  # (SUSY scale)
  1  1     5.45778711E-04   # Y_d(Q)^DRbar
  2  2     1.03697626E-02   # Y_s(Q)^DRbar
  3  3     6.32217941E-01   # Y_b(Q)^DRbar
Block Ye Q=  2.44796300E+03  # (SUSY scale)
  1  1     1.72799619E-04   # Y_e(Q)^DRbar
  2  2     3.57293883E-02   # Y_mu(Q)^DRbar
  3  3     5.97651229E-01   # Y_tau(Q)^DRbar
Block Au Q=  2.44796300E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3    -4.43488000E+03   # A_t(Q)^DRbar
Block Ad Q=  2.44796300E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     2.20173000E+03   # A_b(Q)^DRbar
Block Ae Q=  2.44796300E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     3.64185000E+03   # A_tau(Q)^DRbar
Block MSOFT Q=  2.44796300E+03  # soft SUSY breaking masses at Q
   1    1.62919000E+03  # M_1
   2    3.61340000E+03  # M_2
   3    1.72078000E+03  # M_3
  21    5.13334051E+06  # M^2_(H,d)
  22   -6.88100785E+06  # M^2_(H,u)
  31    3.04231000E+03  # M_(L,11)
  32    3.04231000E+03  # M_(L,22)
  33    3.32718001E+03  # M_(L,33)
  34    1.81372000E+03  # M_(E,11)
  35    1.81372000E+03  # M_(E,22)
  36    3.71686001E+03  # M_(E,33)
  41    2.22692000E+03  # M_(Q,11)
  42    2.22692000E+03  # M_(Q,22)
  43    1.51321000E+03  # M_(Q,33)
  44    3.65349999E+03  # M_(U,11)
  45    3.65349999E+03  # M_(U,22)
  46    3.96489000E+03  # M_(U,33)
  47    3.00973000E+03  # M_(D,11)
  48    3.00973000E+03  # M_(D,22)
  49    3.87152999E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73200000E+02  # m_t(pole)
        23     9.11876000E+01  # m_Z(pole)
        24     8.04071130E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.27589884E+02  # h0
        35     3.45124750E+03  # H0
        36     3.45134000E+03  # A0
        37     3.45242441E+03  # H+
   1000001     2.31470652E+03  # ~d_L
   2000001     3.05938658E+03  # ~d_R
   1000002     2.31341244E+03  # ~u_L
   2000002     3.70204510E+03  # ~u_R
   1000003     2.31470432E+03  # ~s_L
   2000003     3.05938153E+03  # ~s_R
   1000004     2.31341152E+03  # ~c_L
   2000004     3.70204485E+03  # ~c_R
   1000005     1.65494023E+03  # ~b_1
   2000005     3.89666633E+03  # ~b_2
   1000006     1.65202023E+03  # ~t_1
   2000006     3.95995959E+03  # ~t_2
   1000011     3.06265268E+03  # ~e_L-
   2000011     1.81845983E+03  # ~e_R-
   1000012     3.06120287E+03  # ~nu_eL
   1000013     3.06262781E+03  # ~mu_L-
   2000013     1.81842960E+03  # ~mu_R-
   1000014     3.06117093E+03  # ~nu_muL
   1000015     3.34079143E+03  # ~tau_1-
   2000015     3.71080804E+03  # ~tau_2-
   1000016     3.34354864E+03  # ~nu_tauL
   1000021     1.94486809E+03  # ~g
   1000022     1.62436511E+03  # ~chi_10
   1000023     2.60677656E+03  # ~chi_20
   1000025    -2.61029695E+03  # ~chi_30
   1000035     3.56969231E+03  # ~chi_40
   1000024     2.60658647E+03  # ~chi_1+
   1000037     3.56975540E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -1.88290510E-02   # alpha
Block Hmix Q=  2.44796300E+03  # Higgs mixing parameters
   1    2.61534000E+03  # mu
   2    5.31300000E+01  # tan[beta](Q)
   3    2.43793092E+02  # v(Q)
   4    1.22772160E+07  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     9.98742760E-01   # Re[R_st(1,1)]
   1  2     5.01288204E-02   # Re[R_st(1,2)]
   2  1    -5.01288204E-02   # Re[R_st(2,1)]
   2  2     9.98742760E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.99759232E-01   # Re[R_sb(1,1)]
   1  2     2.19426132E-02   # Re[R_sb(1,2)]
   2  1    -2.19426132E-02   # Re[R_sb(2,1)]
   2  2     9.99759232E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     9.94910686E-01   # Re[R_sta(1,1)]
   1  2     1.00760743E-01   # Re[R_sta(1,2)]
   2  1    -1.00760743E-01   # Re[R_sta(2,1)]
   2  2     9.94910686E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.99445555E-01   # Re[N(1,1)]
   1  2     8.02199699E-04   # Re[N(1,2)]
   1  3    -2.80209857E-02   # Re[N(1,3)]
   1  4     1.79656162E-02   # Re[N(1,4)]
   2  1    -3.25113784E-02   # Re[N(2,1)]
   2  2    -5.99167191E-02   # Re[N(2,2)]
   2  3     7.05809459E-01   # Re[N(2,3)]
   2  4    -7.05114178E-01   # Re[N(2,4)]
   3  1    -7.10044963E-03   # Re[N(3,1)]
   3  2     8.95500984E-03   # Re[N(3,2)]
   3  3     7.06928835E-01   # Re[N(3,3)]
   3  4     7.07192346E-01   # Re[N(3,4)]
   4  1     1.08462819E-03   # Re[N(4,1)]
   4  2    -9.98162888E-01   # Re[N(4,2)]
   4  3    -3.60479349E-02   # Re[N(4,3)]
   4  4     4.86848941E-02   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -5.13564123E-02   # Re[U(1,1)]
   1  2     9.98680389E-01   # Re[U(1,2)]
   2  1     9.98680389E-01   # Re[U(2,1)]
   2  2     5.13564123E-02   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -6.91058767E-02   # Re[V(1,1)]
   1  2     9.97609331E-01   # Re[V(1,2)]
   2  1     9.97609331E-01   # Re[V(2,1)]
   2  2     6.91058767E-02   # Re[V(2,2)]
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   2000011     3.90623055E-01   # ~e^-_R
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000011     2.08524460E+00   # ~e^-_L
#    BR                NDA      ID1      ID2
     9.94877425E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     2.75512698E-03    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     2.35651806E-03    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     3.90515111E-01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
DECAY   1000013     2.12706347E+00   # ~mu^-_L
#    BR                NDA      ID1      ID2
     9.75312769E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     4.07736820E-03    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     1.37964618E-03    2     1000025        13   # BR(~mu^-_L -> chi^0_3 mu^-)
     2.29641212E-03    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     8.35304866E-03    2     2000013        23   # BR(~mu^-_L -> ~mu^-_R Z)
     8.58075582E-03    2     2000013        25   # BR(~mu^-_L -> ~mu^-_R h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     6.39719183E+00   # ~tau^-_1
#    BR                NDA      ID1      ID2
     4.18599169E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     2.87074837E-01    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     2.80703152E-01    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.36228416E-02    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
DECAY   2000015     3.00249031E+01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     4.20636183E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     1.11483999E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     1.12298312E-01    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     2.20056101E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     1.21795612E-04    2    -1000037        16   # BR(~tau^-_2 -> chi^-_2 nu_tau)
     6.61774614E-02    2     1000016       -24   # BR(~tau^-_2 -> ~nu_tau W^-)
     3.38723320E-02    2     1000015        23   # BR(~tau^-_2 -> ~tau^-_1 Z)
     3.52930897E-02    2     1000015        25   # BR(~tau^-_2 -> ~tau^-_1 h^0)
DECAY   1000012     2.09425761E+00   # ~nu_e
#    BR                NDA      ID1      ID2
     9.94947218E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     7.52984808E-04    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     4.22546225E-03    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
DECAY   1000014     2.13441289E+00   # ~nu_mu
#    BR                NDA      ID1      ID2
     9.76202852E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     7.38729628E-04    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     6.89024092E-03    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     1.60956267E-02    2     2000013        24   # BR(~nu_mu -> ~mu^-_R W^+)
DECAY   1000016     6.24180698E+00   # ~nu_tau
#    BR                NDA      ID1      ID2
     4.12311867E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     5.61724838E-04    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     5.87070914E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000001     1.11880476E+01   # ~d_L
#    BR                NDA      ID1      ID2
     7.85275664E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     9.92147243E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000001     6.11798679E+01   # ~d_R
#    BR                NDA      ID1      ID2
     1.50716948E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.84925657E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000003     1.11879439E+01   # ~s_L
#    BR                NDA      ID1      ID2
     7.85308786E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     9.92146912E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   2000003     6.11842592E+01   # ~s_R
#    BR                NDA      ID1      ID2
     1.50704863E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     9.84848783E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000005     3.24639908E-04   # ~b_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
DECAY   2000005     1.62405865E+02   # ~b_2
#    BR                NDA      ID1      ID2
     9.60493178E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     2.89668241E-02    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     2.89877890E-02    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     5.75879290E-02    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     7.50369389E-01    2     1000021         5   # BR(~b_2 -> ~g b)
     5.96555700E-02    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
     3.07904301E-02    2     1000005        23   # BR(~b_2 -> ~b_1 Z)
     3.40011634E-02    2     1000005        25   # BR(~b_2 -> ~b_1 h^0)
DECAY   1000002     1.11206044E+01   # ~u_L
#    BR                NDA      ID1      ID2
     7.74760275E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     9.92252397E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000002     1.13289706E+02   # ~u_R
#    BR                NDA      ID1      ID2
     4.98127067E-02    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     9.50165768E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000004     1.11204813E+01   # ~c_L
#    BR                NDA      ID1      ID2
     7.74767886E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     9.92252321E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   2000004     1.13290090E+02   # ~c_R
#    BR                NDA      ID1      ID2
     4.98125077E-02    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     9.50161942E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000006     7.17933023E-11   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     2.76225827E+02   # ~t_2
#    BR                NDA      ID1      ID2
     2.28607310E-02    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     3.13631300E-02    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     3.07445964E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     6.24568813E-02    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     4.52231007E-01    2     1000021         6   # BR(~t_2 -> ~g t)
     1.96137072E-01    2     1000005        24   # BR(~t_2 -> ~b_1 W^+)
     1.01287981E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     1.02802964E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     3.20663092E+01   # chi^+_1
#    BR                NDA      ID1      ID2
     2.72592655E-04    2    -2000013        14   # BR(chi^+_1 -> ~mu^+_R nu_mu)
     2.10150567E-04    2     1000002        -1   # BR(chi^+_1 -> ~u_L d_bar)
     2.21879333E-04    2     1000004        -3   # BR(chi^+_1 -> ~c_L s_bar)
     3.52118346E-01    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     1.15082554E-04    2    -1000001         2   # BR(chi^+_1 -> ~d^*_L u)
     1.16590100E-04    2    -1000003         4   # BR(chi^+_1 -> ~s^*_L c)
     6.02948329E-01    2    -1000005         6   # BR(chi^+_1 -> ~b^*_1 t)
     4.30693852E-02    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     8.63065830E-04    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
DECAY   1000037     1.35078359E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     7.34992019E-03    2    -1000011        12   # BR(chi^+_2 -> ~e^+_L nu_e)
     7.35064753E-03    2    -1000013        14   # BR(chi^+_2 -> ~mu^+_L nu_mu)
     1.59442209E-03    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     7.37293586E-03    2     1000012       -11   # BR(chi^+_2 -> ~nu_e e^+)
     7.37384346E-03    2     1000014       -13   # BR(chi^+_2 -> ~nu_mu mu^+)
     1.58805455E-03    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     1.06261610E-01    2     1000002        -1   # BR(chi^+_2 -> ~u_L d_bar)
     1.06261790E-01    2     1000004        -3   # BR(chi^+_2 -> ~c_L s_bar)
     1.93248209E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     1.06317393E-01    2    -1000001         2   # BR(chi^+_2 -> ~d^*_L u)
     1.06317217E-01    2    -1000003         4   # BR(chi^+_2 -> ~s^*_L c)
     1.91495066E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     3.80256516E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     3.85801057E-02    2     1000025        24   # BR(chi^+_2 -> chi^0_3 W^+)
     3.99505342E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     3.82817324E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     8.82296028E-04    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     7.57725426E-04    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
     7.51422726E-04    3     1000025        -5         6   # BR(chi^+_2 -> chi^0_3 b_bar t)
     2.22354385E-04    3     1000021        -5         6   # BR(chi^+_2 -> ~g b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     3.24834714E+01   # chi^0_2
#    BR                NDA      ID1      ID2
     1.93165132E-04    2     2000013       -13   # BR(chi^0_2 -> ~mu^-_R mu^+)
     1.93165132E-04    2    -2000013        13   # BR(chi^0_2 -> ~mu^+_R mu^-)
     3.05451629E-01    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     3.05451629E-01    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     1.70863921E-01    2     1000005        -5   # BR(chi^0_2 -> ~b_1 b_bar)
     1.70863921E-01    2    -1000005         5   # BR(chi^0_2 -> ~b^*_1 b)
     2.04265974E-03    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     4.37254019E-02    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.99617896E-04    3     1000022         6        -6   # BR(chi^0_2 -> chi^0_1 t t_bar)
DECAY   1000025     3.12463994E+01   # chi^0_3
#    BR                NDA      ID1      ID2
     1.43211442E-04    2     2000013       -13   # BR(chi^0_3 -> ~mu^-_R mu^+)
     1.43211442E-04    2    -2000013        13   # BR(chi^0_3 -> ~mu^+_R mu^-)
     2.98545866E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     2.98545866E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     1.77989781E-01    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     1.77989781E-01    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     4.37475322E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     2.26134528E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     5.52441600E-04    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
DECAY   1000035     1.35260013E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     3.66086772E-03    2     1000011       -11   # BR(chi^0_4 -> ~e^-_L e^+)
     3.66086772E-03    2    -1000011        11   # BR(chi^0_4 -> ~e^+_L e^-)
     3.66125530E-03    2     1000013       -13   # BR(chi^0_4 -> ~mu^-_L mu^+)
     3.66125530E-03    2    -1000013        13   # BR(chi^0_4 -> ~mu^+_L mu^-)
     7.95174985E-04    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     7.95174985E-04    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     3.68940855E-03    2     1000012       -12   # BR(chi^0_4 -> ~nu_e nu_bar_e)
     3.68940855E-03    2    -1000012        12   # BR(chi^0_4 -> ~nu^*_e nu_e)
     3.68983660E-03    2     1000014       -14   # BR(chi^0_4 -> ~nu_mu nu_bar_mu)
     3.68983660E-03    2    -1000014        14   # BR(chi^0_4 -> ~nu^*_mu nu_mu)
     7.93196568E-04    2     1000016       -16   # BR(chi^0_4 -> ~nu_tau nu_bar_tau)
     7.93196568E-04    2    -1000016        16   # BR(chi^0_4 -> ~nu^*_tau nu_tau)
     5.30926290E-02    2     1000002        -2   # BR(chi^0_4 -> ~u_L u_bar)
     5.30926290E-02    2    -1000002         2   # BR(chi^0_4 -> ~u^*_L u)
     5.30926311E-02    2     1000004        -4   # BR(chi^0_4 -> ~c_L c_bar)
     5.30926311E-02    2    -1000004         4   # BR(chi^0_4 -> ~c^*_L c)
     9.48193103E-02    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     9.48193103E-02    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     5.30506761E-02    2     1000001        -1   # BR(chi^0_4 -> ~d_L d_bar)
     5.30506761E-02    2    -1000001         1   # BR(chi^0_4 -> ~d^*_L d)
     5.30506766E-02    2     1000003        -3   # BR(chi^0_4 -> ~s_L s_bar)
     5.30506766E-02    2    -1000003         3   # BR(chi^0_4 -> ~s^*_L s)
     9.73103748E-02    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     9.73103748E-02    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     3.91294705E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     3.91294705E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     8.74750289E-04    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     3.83835807E-02    2     1000025        23   # BR(chi^0_4 -> chi^0_3 Z)
     3.74152559E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
     8.28171400E-04    2     1000025        25   # BR(chi^0_4 -> chi^0_3 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.89680821E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     5.38287802E-04    3     1000025         6        -6   # BR(chi^0_4 -> chi^0_3 t t_bar)
     7.82048571E-04    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     7.82048571E-04    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
     1.43022374E-04    3     1000021         6        -6   # BR(chi^0_4 -> ~g t t_bar)
     1.69841514E-04    3     1000021         5        -5   # BR(chi^0_4 -> ~g b b_bar)
DECAY   1000021     5.51470625E+00   # ~g
#    BR                NDA      ID1      ID2
     2.20954735E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     2.20954735E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     2.79044632E-01    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     2.79044632E-01    2    -1000005         5   # BR(~g -> ~b^*_1 b)
#    BR                NDA      ID1      ID2       ID3
DECAY        25     3.14496503E-03   # h^0
#    BR                NDA      ID1      ID2
     3.63936175E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     9.78389144E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     1.26770884E-04    2           3        -3   # BR(h^0 -> s s_bar)
     4.36430806E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.96153835E-02    2           4        -4   # BR(h^0 -> c c_bar)
     1.05372863E-01    2          21        21   # BR(h^0 -> g g)
     3.82229573E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.41342346E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.41342346E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.41342346E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     4.94698212E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     4.94698212E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.41342346E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.41342346E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.41342346E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     4.94698212E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     4.94698212E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     1.13379252E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     1.13379252E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     1.14054128E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     6.74876499E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     5.26403669E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     5.26403669E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     5.09531757E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     3.91428369E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     4.04925899E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     1.19727605E+02   # H^0
#    BR                NDA      ID1      ID2
     7.31828188E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     2.04763676E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     2.13274833E-04    2           3        -3   # BR(H^0 -> s s_bar)
     7.92750875E-01    2           5        -5   # BR(H^0 -> b b_bar)
     4.82139503E-04    2           6        -6   # BR(H^0 -> t t_bar)
     9.54590101E-04    2     1000006  -1000006   # BR(H^0 -> ~t_1 ~t^*_1)
DECAY        36     1.19607214E+02   # A^0
#    BR                NDA      ID1      ID2
     7.32584707E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     2.04975565E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     2.13495374E-04    2           3        -3   # BR(A^0 -> s s_bar)
     7.93569755E-01    2           5        -5   # BR(A^0 -> b b_bar)
     4.84469717E-04    2           6        -6   # BR(A^0 -> t t_bar)
DECAY        37     1.06694838E+02   # H^+
#    BR                NDA      ID1      ID2
     8.21501290E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     2.29854226E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.07602802E-04    2          -3         4   # BR(H^+ -> s_bar c)
     7.68218099E-01    2          -5         6   # BR(H^+ -> b_bar t)
     8.97966088E-04    2    -1000005   1000006   # BR(H^+ -> ~b^*_1 ~t_1)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00101134E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    2.82279589E+03    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    2.82279690E+03        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99999642E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    3.54616849E-04    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    3.54258572E-04        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00101134E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    2.82279589E+03    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    2.82279690E+03        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    1.00000000E+00        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    9.05392292E-11        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    1.00000000E+00        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    9.05392292E-11        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.16228168E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    1.88784040E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.47194698E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    9.05392292E-11        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    1.00000000E+00        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.09583464E-04   # BR(b -> s gamma)
    2    1.59295250E-06   # BR(b -> s mu+ mu-)
    3    3.52982026E-05   # BR(b -> s nu nu)
    4    2.27810489E-15   # BR(Bd -> e+ e-)
    5    9.73140935E-11   # BR(Bd -> mu+ mu-)
    6    2.01099632E-08   # BR(Bd -> tau+ tau-)
    7    7.61461025E-14   # BR(Bs -> e+ e-)
    8    3.25283146E-09   # BR(Bs -> mu+ mu-)
    9    6.81382843E-07   # BR(Bs -> tau+ tau-)
   10    9.62366704E-05   # BR(B_u -> tau nu)
   11    9.94102792E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41298633E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93381166E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15575953E-03   # epsilon_K
   17    2.28160792E-15   # Delta(M_K)
   18    2.48292017E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.29765219E-11   # BR(K^+ -> pi^+ nu nu)
   20    1.60759863E-15   # Delta(g-2)_electron/2
   21    6.87247094E-11   # Delta(g-2)_muon/2
   22    1.84474944E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    6.87140483E-05   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88852617E-01   # C7
     0305 4422   00   2    -1.75335344E-01   # C7
     0305 4322   00   2     1.78553174E-04   # C7'
     0305 6421   00   0    -9.52384852E-02   # C8
     0305 6421   00   2    -9.16353812E-02   # C8
     0305 6321   00   2     4.13223346E-06   # C8'
 03051111 4133   00   0     1.16465527E+00   # C9 e+e-
 03051111 4133   00   2     1.16496080E+00   # C9 e+e-
 03051111 4233   00   2     5.07332922E-05   # C9' e+e-
 03051111 4137   00   0    -3.98903975E+00   # C10 e+e-
 03051111 4137   00   2    -3.98962989E+00   # C10 e+e-
 03051111 4237   00   2    -1.45687341E-03   # C10' e+e-
 03051313 4133   00   0     1.16465527E+00   # C9 mu+mu-
 03051313 4133   00   2     1.16495976E+00   # C9 mu+mu-
 03051313 4233   00   2     5.07320005E-05   # C9' mu+mu-
 03051313 4137   00   0    -3.98903975E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.98963093E+00   # C10 mu+mu-
 03051313 4237   00   2    -1.45687476E-03   # C10' mu+mu-
 03051212 4137   00   0     1.50583348E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50597560E+00   # C11 nu_1 nu_1
 03051212 4237   00   2     3.51541401E-04   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50583348E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50597560E+00   # C11 nu_2 nu_2
 03051414 4237   00   2     3.51541521E-04   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50583348E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50597570E+00   # C11 nu_3 nu_3
 03051616 4237   00   2     3.51575438E-04   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.86112808E-07   # C7
     0305 4422   00   2     2.80739205E-05   # C7
     0305 4322   00   2     4.13365151E-06   # C7'
     0305 6421   00   0     3.30731046E-07   # C8
     0305 6421   00   2    -3.63231283E-05   # C8
     0305 6321   00   2     1.31831453E-06   # C8'
 03051111 4133   00   2     9.45700883E-07   # C9 e+e-
 03051111 4233   00   2     1.14447130E-06   # C9' e+e-
 03051111 4137   00   2    -2.82754377E-06   # C10 e+e-
 03051111 4237   00   2    -3.28690318E-05   # C10' e+e-
 03051313 4133   00   2     9.45685804E-07   # C9 mu+mu-
 03051313 4233   00   2     1.14445742E-06   # C9' mu+mu-
 03051313 4137   00   2    -2.82755878E-06   # C10 mu+mu-
 03051313 4237   00   2    -3.28690651E-05   # C10' mu+mu-
 03051212 4137   00   2     6.77537057E-07   # C11 nu_1 nu_1
 03051212 4237   00   2     7.93124811E-06   # C11' nu_1 nu_1
 03051414 4137   00   2     6.77536013E-07   # C11 nu_2 nu_2
 03051414 4237   00   2     7.93124811E-06   # C11' nu_2 nu_2
 03051616 4137   00   2     6.79217561E-07   # C11 nu_3 nu_3
 03051616 4237   00   2     7.93124807E-06   # C11' nu_3 nu_3
