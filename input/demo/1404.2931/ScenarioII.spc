# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 09.06.2023,  11:30
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
    3    6.84256288E+00  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    7.05400000E+02  # M_1
    2    1.92840000E+03  # M_2
    3    1.42700000E+03  # M_3
   23    8.42300000E+02  # mu
   25    6.60000000E+00  # tan(beta)
   26    1.56690000E+03  # m_A, pole mass
   31    2.93320000E+03  # M^2_L11
   32    2.93320000E+03  # M^2_L22
   33    2.93320000E+03  # M^2_L33
   34    2.93320000E+03  # M^2_E11
   35    2.93320000E+03  # M^2_E22
   36    2.93320000E+03  # M^2_E33
   41    1.23850000E+03  # M^2_Q11
   42    1.23850000E+03  # M^2_Q22
   43    2.35210000E+03  # M^2_Q33
   44    1.23850000E+03  # M^2_U11
   45    1.23850000E+03  # M^2_U22
   46    7.74100000E+02  # M^2_U33
   47    1.23850000E+03  # M^2_D11
   48    1.23850000E+03  # M^2_D22
   49    2.35210000E+03  # M^2_D33
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
   1    3.61705932E-01  # g'(Q)^DRbar
   2    6.34792431E-01  # g(Q)^DRbar
   3    1.04970860E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.78654250E-06   # Y_u(Q)^DRbar
  2  2     3.71963581E-03   # Y_c(Q)^DRbar
  3  3     8.57075733E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     9.19905931E-05   # Y_d(Q)^DRbar
  2  2     1.93180146E-03   # Y_s(Q)^DRbar
  3  3     9.70307739E-02   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.98177992E-05   # Y_e(Q)^DRbar
  2  2     4.09767696E-03   # Y_mu(Q)^DRbar
  3  3     6.89113360E-02   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3    -3.70399006E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     0.00000000E+00   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     0.00000000E+00   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    7.05400000E+02  # M_1
   2    1.92840000E+03  # M_2
   3    1.42700000E+03  # M_3
  21    1.67547947E+06  # M^2_(H,d)
  22   -8.96009382E+05  # M^2_(H,u)
  31    2.93320000E+03  # M_(L,11)
  32    2.93320000E+03  # M_(L,22)
  33    2.93320000E+03  # M_(L,33)
  34    2.93320000E+03  # M_(E,11)
  35    2.93320000E+03  # M_(E,22)
  36    2.93320000E+03  # M_(E,33)
  41    1.23850000E+03  # M_(Q,11)
  42    1.23850000E+03  # M_(Q,22)
  43    2.35210000E+03  # M_(Q,33)
  44    1.23850000E+03  # M_(U,11)
  45    1.23850000E+03  # M_(U,22)
  46    7.74100000E+02  # M_(U,33)
  47    1.23850000E+03  # M_(D,11)
  48    1.23850000E+03  # M_(D,22)
  49    2.35210000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02910239E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.25195498E+02  # h0
        35     1.56706614E+03  # H0
        36     1.56690000E+03  # A0
        37     1.56898009E+03  # H+
   1000001     1.28614428E+03  # ~d_L
   2000001     1.28142601E+03  # ~d_R
   1000002     1.28400308E+03  # ~u_L
   2000002     1.27800638E+03  # ~u_R
   1000003     1.28614682E+03  # ~s_L
   2000003     1.28142334E+03  # ~s_R
   1000004     1.28400362E+03  # ~c_L
   2000004     1.27800607E+03  # ~c_R
   1000005     2.36191732E+03  # ~b_1
   2000005     2.37428262E+03  # ~b_2
   1000006     8.74812179E+02  # ~t_1
   2000006     2.38036373E+03  # ~t_2
   1000011     2.93500091E+03  # ~e_L-
   2000011     2.93695812E+03  # ~e_R-
   1000012     2.93366810E+03  # ~nu_eL
   1000013     2.93499593E+03  # ~mu_L-
   2000013     2.93696363E+03  # ~mu_R-
   1000014     2.93366828E+03  # ~nu_muL
   1000015     2.93408839E+03  # ~tau_1-
   2000015     2.93802219E+03  # ~tau_2-
   1000016     2.93371847E+03  # ~nu_tauL
   1000021     1.49676436E+03  # ~g
   1000022     6.98911340E+02  # ~chi_10
   1000023    -8.50532972E+02  # ~chi_20
   1000025     8.53956363E+02  # ~chi_30
   1000035     1.94020143E+03  # ~chi_40
   1000024     8.45589627E+02  # ~chi_1+
   1000037     1.94041876E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -1.51895004E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    8.42300000E+02  # mu
   2    6.60000000E+00  # tan[beta](Q)
   3    2.44355062E+02  # v(Q)
   4    2.48024920E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     1.18489092E-01   # Re[R_st(1,1)]
   1  2     9.92955354E-01   # Re[R_st(1,2)]
   2  1    -9.92955354E-01   # Re[R_st(2,1)]
   2  2     1.18489092E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.74524103E-01   # Re[R_sb(1,1)]
   1  2     2.24282794E-01   # Re[R_sb(1,2)]
   2  1    -2.24282794E-01   # Re[R_sb(2,1)]
   2  2     9.74524103E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     8.70712172E-01   # Re[R_sta(1,1)]
   1  2     4.91792958E-01   # Re[R_sta(1,2)]
   2  1    -4.91792958E-01   # Re[R_sta(2,1)]
   2  2     8.70712172E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.68975327E-01   # Re[N(1,1)]
   1  2     1.18218029E-02   # Re[N(1,2)]
   1  3    -1.86544045E-01   # Re[N(1,3)]
   1  4     1.61704608E-01   # Re[N(1,4)]
   2  1    -1.75876027E-02   # Re[N(2,1)]
   2  2     1.65858265E-02   # Re[N(2,2)]
   2  3     7.05977892E-01   # Re[N(2,3)]
   2  4     7.07821165E-01   # Re[N(2,4)]
   3  1    -2.46522308E-01   # Re[N(3,1)]
   3  2    -5.59911743E-02   # Re[N(3,2)]
   3  3     6.82619182E-01   # Re[N(3,3)]
   3  4    -6.85655010E-01   # Re[N(3,4)]
   4  1    -2.05999285E-03   # Re[N(4,1)]
   4  2     9.98223494E-01   # Re[N(4,2)]
   4  3     2.87678159E-02   # Re[N(4,3)]
   4  4    -5.21346858E-02   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -4.07236348E-02   # Re[U(1,1)]
   1  2     9.99170449E-01   # Re[U(1,2)]
   2  1     9.99170449E-01   # Re[U(2,1)]
   2  2     4.07236348E-02   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -7.38101240E-02   # Re[V(1,1)]
   1  2     9.97272313E-01   # Re[V(1,2)]
   2  1     9.97272313E-01   # Re[V(2,1)]
   2  2     7.38101240E-02   # Re[V(2,2)]
DECAY   1000011     1.46197296E+01   # ~e^-_L
#    BR                NDA      ID1      ID2
     2.09026021E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     2.60225755E-02    2     1000025        11   # BR(~e^-_L -> chi^0_3 e^-)
     2.53565142E-01    2     1000035        11   # BR(~e^-_L -> chi^0_4 e^-)
     2.24434660E-03    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
     5.09112818E-01    2    -1000037        12   # BR(~e^-_L -> chi^-_2 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000011     1.35575958E+01   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.42270858E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
     2.92764948E-04    2     1000023        11   # BR(~e^-_R -> chi^0_2 e^-)
     5.74348057E-02    2     1000025        11   # BR(~e^-_R -> chi^0_3 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000013     1.46182661E+01   # ~mu^-_L
#    BR                NDA      ID1      ID2
     2.10831494E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     2.60854654E-02    2     1000025        13   # BR(~mu^-_L -> chi^0_3 mu^-)
     2.52928162E-01    2     1000035        13   # BR(~mu^-_L -> chi^0_4 mu^-)
     2.27497899E-03    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     5.07833874E-01    2    -1000037        14   # BR(~mu^-_L -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.35615354E+01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.40078105E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     3.34934753E-04    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     5.74096791E-02    2     1000025        13   # BR(~mu^-_R -> chi^0_3 mu^-)
     7.15619465E-04    2     1000035        13   # BR(~mu^-_R -> chi^0_4 mu^-)
     1.43372931E-03    2    -1000037        14   # BR(~mu^-_R -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     1.47161273E+01   # ~tau^-_1
#    BR                NDA      ID1      ID2
     3.77605163E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     6.30966416E-03    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     3.46537924E-02    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.89784711E-01    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     1.05975106E-02    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
     3.81049158E-01    2    -1000037        16   # BR(~tau^-_1 -> chi^-_2 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000015     1.41616393E+01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     7.26711564E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     1.01400978E-02    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     6.11523650E-02    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     6.45859523E-02    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     7.74756368E-03    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     1.29662457E-01    2    -1000037        16   # BR(~tau^-_2 -> chi^-_2 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     1.46309419E+01   # ~nu_e
#    BR                NDA      ID1      ID2
     2.27417496E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     4.77356153E-04    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     4.80488203E-03    2     1000025        12   # BR(~nu_e -> chi^0_3 nu_e)
     2.54090914E-01    2     1000035        12   # BR(~nu_e -> chi^0_4 nu_e)
     7.36252870E-03    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
     5.05846823E-01    2     1000037        11   # BR(~nu_e -> chi^+_2 e^-)
DECAY   1000014     1.46317679E+01   # ~nu_mu
#    BR                NDA      ID1      ID2
     2.27404676E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     4.77329248E-04    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     4.80461121E-03    2     1000025        14   # BR(~nu_mu -> chi^0_3 nu_mu)
     2.54076635E-01    2     1000035        14   # BR(~nu_mu -> chi^0_4 nu_mu)
     7.41829947E-03    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     5.05818449E-01    2     1000037        13   # BR(~nu_mu -> chi^+_2 mu^-)
DECAY   1000016     1.48645375E+01   # ~nu_tau
#    BR                NDA      ID1      ID2
     2.23848408E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     4.69865550E-04    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     4.72948479E-03    2     1000025        16   # BR(~nu_tau -> chi^0_3 nu_tau)
     2.50115535E-01    2     1000035        16   # BR(~nu_tau -> chi^0_4 nu_tau)
     2.28891804E-02    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     4.97947526E-01    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
DECAY   2000001     3.57442954E-01   # ~d_R
#    BR                NDA      ID1      ID2
     9.60857156E-01    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     2.00804920E-04    2     1000023         1   # BR(~d_R -> chi^0_2 d)
     3.89417050E-02    2     1000025         1   # BR(~d_R -> chi^0_3 d)
DECAY   1000001     1.04507868E-01   # ~d_L
#    BR                NDA      ID1      ID2
     9.39762298E-01    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     6.20165157E-03    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.29725716E-03    2     1000025         1   # BR(~d_L -> chi^0_3 d)
     5.27387935E-02    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
DECAY   2000003     3.57405192E-01   # ~s_R
#    BR                NDA      ID1      ID2
     9.60716692E-01    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     2.25838633E-04    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     3.89100030E-02    2     1000025         3   # BR(~s_R -> chi^0_3 s)
     1.47466487E-04    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
DECAY   1000003     1.04745666E-01   # ~s_L
#    BR                NDA      ID1      ID2
     9.38477522E-01    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     6.38752906E-03    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.68375714E-03    2     1000025         3   # BR(~s_L -> chi^0_3 s)
     5.34511917E-02    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
DECAY   1000005     1.47152922E+02   # ~b_1
#    BR                NDA      ID1      ID2
     2.93301145E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     1.06727721E-03    2     1000023         5   # BR(~b_1 -> chi^0_2 b)
     8.18229119E-04    2     1000025         5   # BR(~b_1 -> chi^0_3 b)
     6.42655065E-03    2     1000035         5   # BR(~b_1 -> chi^0_4 b)
     1.65412384E-01    2    -1000024         6   # BR(~b_1 -> chi^-_1 t)
     1.24964553E-02    2    -1000037         6   # BR(~b_1 -> chi^-_2 t)
     3.36549431E-01    2     1000021         5   # BR(~b_1 -> ~g b)
     4.74296661E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     5.73182799E+01   # ~b_2
#    BR                NDA      ID1      ID2
     1.69377738E-02    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     3.16764264E-03    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     4.48148053E-03    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     9.69918986E-04    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     2.89019588E-02    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     1.88897863E-03    2    -1000037         6   # BR(~b_2 -> chi^-_2 t)
     8.77807619E-01    2     1000021         5   # BR(~b_2 -> ~g b)
     6.58446277E-02    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   2000002     1.41926888E+00   # ~u_R
#    BR                NDA      ID1      ID2
     9.61009003E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     1.99995634E-04    2     1000023         2   # BR(~u_R -> chi^0_2 u)
     3.87910008E-02    2     1000025         2   # BR(~u_R -> chi^0_3 u)
DECAY   1000002     1.10782357E-01   # ~u_L
#    BR                NDA      ID1      ID2
     6.82363018E-01    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     2.56712501E-03    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.52732228E-01    2     1000025         2   # BR(~u_L -> chi^0_3 u)
     1.62337630E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
DECAY   2000004     1.41944741E+00   # ~c_R
#    BR                NDA      ID1      ID2
     9.60845614E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     2.41866558E-04    2     1000023         4   # BR(~c_R -> chi^0_2 c)
     3.88171945E-02    2     1000025         4   # BR(~c_R -> chi^0_3 c)
DECAY   1000004     1.10958816E-01   # ~c_L
#    BR                NDA      ID1      ID2
     6.81797993E-01    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     3.03250614E-03    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.53044961E-01    2     1000025         4   # BR(~c_L -> chi^0_3 c)
     1.62124540E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
DECAY   1000006     7.07199054E-02   # ~t_1
#    BR                NDA      ID1      ID2
     2.36186019E-01    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     7.63813926E-01    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000006     1.51648007E+02   # ~t_2
#    BR                NDA      ID1      ID2
     4.95366863E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     8.47316096E-02    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     8.28670642E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     6.96045370E-03    2     1000035         6   # BR(~t_2 -> chi^0_4 t)
     2.92807293E-03    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     1.42223063E-02    2     1000037         5   # BR(~t_2 -> chi^+_2 b)
     3.13996567E-01    2     1000021         6   # BR(~t_2 -> ~g t)
     2.52649058E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     2.36691199E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     2.17361617E-01   # chi^+_1
#    BR                NDA      ID1      ID2
     9.99999996E-01    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000037     4.69338664E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     1.56243254E-01    2     1000002        -1   # BR(chi^+_2 -> ~u_L d_bar)
     1.56233146E-01    2     1000004        -3   # BR(chi^+_2 -> ~c_L s_bar)
     1.16609183E-04    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     1.56023864E-01    2    -1000001         2   # BR(chi^+_2 -> ~d^*_L u)
     1.55934981E-01    2    -1000003         4   # BR(chi^+_2 -> ~s^*_L c)
     3.34847690E-03    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     9.02045581E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     8.67064421E-02    2     1000025        24   # BR(chi^+_2 -> chi^0_3 W^+)
     9.27591547E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     9.51857296E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     2.70742176E-03    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     2.18286137E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
     2.10279222E-03    3     1000025        -5         6   # BR(chi^+_2 -> chi^0_3 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     2.09335265E-01   # chi^0_2
#    BR                NDA      ID1      ID2
     9.98268578E-01    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     1.73132155E-03    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000025     1.99401169E-01   # chi^0_3
#    BR                NDA      ID1      ID2
     1.32418956E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     9.86757283E-01    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000035     4.69144941E+01   # chi^0_4
#    BR                NDA      ID1      ID2
     7.82056038E-02    2     1000002        -2   # BR(chi^0_4 -> ~u_L u_bar)
     7.82056038E-02    2    -1000002         2   # BR(chi^0_4 -> ~u^*_L u)
     7.82003066E-02    2     1000004        -4   # BR(chi^0_4 -> ~c_L c_bar)
     7.82003066E-02    2    -1000004         4   # BR(chi^0_4 -> ~c^*_L c)
     7.79212641E-02    2     1000001        -1   # BR(chi^0_4 -> ~d_L d_bar)
     7.79212641E-02    2    -1000001         1   # BR(chi^0_4 -> ~d^*_L d)
     7.78771062E-02    2     1000003        -3   # BR(chi^0_4 -> ~s_L s_bar)
     7.78771062E-02    2    -1000003         3   # BR(chi^0_4 -> ~s^*_L s)
     9.08453436E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     9.08453436E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     3.09444974E-04    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     8.50035854E-02    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     6.82079572E-03    2     1000025        23   # BR(chi^0_4 -> chi^0_3 Z)
     3.23840629E-03    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     7.22395803E-03    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
     8.46494145E-02    2     1000025        25   # BR(chi^0_4 -> chi^0_3 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.38492921E-03    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     9.71371554E-04    3     1000025         6        -6   # BR(chi^0_4 -> chi^0_3 t t_bar)
     1.97546920E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     1.97546920E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     3.13465124E+01   # ~g
#    BR                NDA      ID1      ID2
     3.84211322E-02    2     2000002        -2   # BR(~g -> ~u_R u_bar)
     3.84211322E-02    2    -2000002         2   # BR(~g -> ~u^*_R u)
     3.65008227E-02    2     1000002        -2   # BR(~g -> ~u_L u_bar)
     3.65008227E-02    2    -1000002         2   # BR(~g -> ~u^*_L u)
     3.84166378E-02    2     2000004        -4   # BR(~g -> ~c_R c_bar)
     3.84166378E-02    2    -2000004         4   # BR(~g -> ~c^*_R c)
     3.65040413E-02    2     1000004        -4   # BR(~g -> ~c_L c_bar)
     3.65040413E-02    2    -1000004         4   # BR(~g -> ~c^*_L c)
     2.03742573E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     2.03742573E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     3.73211394E-02    2     2000001        -1   # BR(~g -> ~d_R d_bar)
     3.73211394E-02    2    -2000001         1   # BR(~g -> ~d^*_R d)
     3.58249503E-02    2     1000001        -1   # BR(~g -> ~d_L d_bar)
     3.58249503E-02    2    -1000001         1   # BR(~g -> ~d^*_L d)
     3.73210546E-02    2     2000003        -3   # BR(~g -> ~s_R s_bar)
     3.73210546E-02    2    -2000003         3   # BR(~g -> ~s^*_R s)
     3.58250623E-02    2     1000003        -3   # BR(~g -> ~s_L s_bar)
     3.58250623E-02    2    -1000003         3   # BR(~g -> ~s^*_L s)
#    BR                NDA      ID1      ID2       ID3
     1.08371715E-04    3     1000023         6        -6   # BR(~g -> chi^0_2 t t_bar)
DECAY        25     3.63413605E-03   # h^0
#    BR                NDA      ID1      ID2
     2.74348377E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     7.74673366E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.38674308E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.88192281E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.56081187E-02    2           4        -4   # BR(h^0 -> c c_bar)
     8.02137424E-02    2          21        21   # BR(h^0 -> g g)
     3.11929852E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     9.67489299E-03    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     9.67489299E-03    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     9.67489299E-03    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     3.38621255E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     3.38621255E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     9.67489299E-03    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     9.67489299E-03    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     9.67489299E-03    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     3.38621255E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     3.38621255E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     7.48565772E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     7.48565772E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     6.41627804E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     4.27751870E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     2.99426309E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     2.99426309E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     2.99426309E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     2.99426309E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     2.99426309E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     3.14783985E+00   # H^0
#    BR                NDA      ID1      ID2
     1.62488046E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     4.59540698E-02    2          15       -15   # BR(H^0 -> tau^- tau^+)
     1.24151222E-04    2           3        -3   # BR(H^0 -> s s_bar)
     3.13218180E-01    2           5        -5   # BR(H^0 -> b b_bar)
     5.50132147E-01    2           6        -6   # BR(H^0 -> t t_bar)
     2.99856036E-03    2     1000022   1000022   # BR(H^0 -> chi^0_1 chi^0_1)
     5.94999488E-02    2     1000022   1000023   # BR(H^0 -> chi^0_1 chi^0_2)
     4.57884635E-04    2     1000022   1000025   # BR(H^0 -> chi^0_1 chi^0_3)
     8.78442181E-03    2          23        23   # BR(H^0 -> Z Z)
     1.67339376E-02    2          24       -24   # BR(H^0 -> W^+ W^-)
     1.60071135E-03    2          25        25   # BR(H^0 -> h^0 h^0)
     3.19899993E-04    2          21        21   # BR(H^0 -> g g)
DECAY        36     3.11478775E+00   # A^0
#    BR                NDA      ID1      ID2
     1.64271071E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     4.64585743E-02    2          15       -15   # BR(A^0 -> tau^- tau^+)
     1.25513409E-04    2           3        -3   # BR(A^0 -> s s_bar)
     3.16654775E-01    2           5        -5   # BR(A^0 -> b b_bar)
     5.62605512E-01    2           6        -6   # BR(A^0 -> t t_bar)
     2.52328528E-02    2     1000022   1000022   # BR(A^0 -> chi^0_1 chi^0_1)
     7.13344929E-04    2     1000022   1000023   # BR(A^0 -> chi^0_1 chi^0_2)
     4.70043769E-02    2     1000022   1000025   # BR(A^0 -> chi^0_1 chi^0_3)
     8.88661265E-04    2          25        23   # BR(A^0 -> h^0 Z)
     1.40198329E-04    2          21        21   # BR(A^0 -> g g)
DECAY        37     2.71304179E+00   # H^+
#    BR                NDA      ID1      ID2
     1.88846619E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     5.34089460E-02    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     1.36632421E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.64643863E-01    2          -5         6   # BR(H^+ -> b_bar t)
     8.06216093E-02    2     1000024   1000022   # BR(H^+ -> chi^+_1 chi^0_1)
     9.99812246E-04    2          25        24   # BR(H^+ -> h^0 W^+)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.02020996E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    4.35397900E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    4.35600000E+01        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99536043E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    2.34207980E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    2.29568411E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.02020996E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    4.35397900E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    4.35600000E+01        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99997679E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    2.32128108E-06        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99997679E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    2.32128108E-06        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.08206544E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    1.92319138E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.67331738E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    2.32128108E-06        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99997679E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    2.98497508E-04   # BR(b -> s gamma)
    2    1.59206521E-06   # BR(b -> s mu+ mu-)
    3    3.54513108E-05   # BR(b -> s nu nu)
    4    5.20077944E-11   # BR(Bd -> mu+ mu-)
    5    3.39740003E-09   # BR(Bs -> mu+ mu-)
    6    1.10499961E-04   # BR(B_u -> tau nu)
    7    9.99190063E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.72968894E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.02594392E+01   # |Delta(M_Bs)| [ps^-1] 
   20    1.18051561E-15   # Delta(g-2)_electron/2
   21    1.78385605E-11   # Delta(g-2)_muon/2
   22   -4.28165008E-09   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    1.46707104E-02   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87790605E-01   # C7
     0305 4422   00   1     1.26458916E-03   # C7
     0305 4322   00   1     8.63168389E-06   # C7'
     0305 6421   00   0    -9.48999873E-02   # C8
     0305 6421   00   1    -6.31035091E-03   # C8
     0305 6321   00   1    -1.38883614E-04   # C8'
 03051111 4133   00   0    -4.71105509E-01   # C9 e+e-
 03051111 4133   00   1     1.30456288E-03   # C9 e+e-
 03051111 4233   00   1    -8.25632119E-07   # C9' e+e-
 03051111 4137   00   0    -3.93833634E+00   # C10 e+e-
 03051111 4137   00   1    -5.35922729E-04   # C10 e+e-
 03051111 4237   00   1     2.33224132E-05   # C10' e+e-
 03051313 4133   00   0    -4.71105509E-01   # C9 mu+mu-
 03051313 4133   00   1     1.30455991E-03   # C9 mu+mu-
 03051313 4233   00   1    -8.25635208E-07   # C9' mu+mu-
 03051313 4137   00   0    -3.93833634E+00   # C10 mu+mu-
 03051313 4137   00   1    -5.35925710E-04   # C10 mu+mu-
 03051313 4237   00   1     2.33224163E-05   # C10' mu+mu-
 03051212 4237   00   0     1.49061874E+00   # C11 nu_1 nu_1
 03051212 4237   00   1     1.25543064E-04   # C11 nu_1 nu_1
 03051212 4137   00   1    -5.62093479E-06   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49061874E+00   # C11 nu_2 nu_2
 03051414 4237   00   1     1.25543064E-04   # C11 nu_2 nu_2
 03051414 4137   00   1    -5.62093419E-06   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49061874E+00   # C11 nu_3 nu_3
 03051616 4237   00   1     1.25543064E-04   # C11 nu_3 nu_3
 03051616 4137   00   1    -5.62076566E-06   # C11' nu_3 nu_3
