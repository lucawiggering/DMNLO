# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 14.03.2023,  13:56
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
    3    6.03028203E+00  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    3.35000000E+02  # M_1
    2    1.95410000E+03  # M_2
    3    1.94560000E+03  # M_3
   11    3.03370000E+03  # A_t
   12    0.00000000E+00  # A_b
   13    0.00000000E+00  # A_tau
   23    2.92580000E+03  # mu
   25    5.80000000E+00  # tan(beta)
   26    9.48800000E+02  # m_A, pole mass
   31    3.26390000E+03  # M^2_L11
   32    3.26390000E+03  # M^2_L22
   33    3.26390000E+03  # M^2_L33
   34    3.26390000E+03  # M^2_E11
   35    3.26390000E+03  # M^2_E22
   36    3.26390000E+03  # M^2_E33
   41    3.21510000E+03  # M^2_Q11
   42    3.21510000E+03  # M^2_Q22
   43    1.57800000E+03  # M^2_Q33
   44    3.21510000E+03  # M^2_U11
   45    3.21510000E+03  # M^2_U22
   46    6.09200000E+02  # M^2_U33
   47    3.21510000E+03  # M^2_D11
   48    3.21510000E+03  # M^2_D22
   49    1.57800000E+03  # M^2_D33
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
   1    3.61096239E-01  # g'(Q)^DRbar
   2    6.33203493E-01  # g(Q)^DRbar
   3    1.03694581E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.70018197E-06   # Y_u(Q)^DRbar
  2  2     3.68307618E-03   # Y_c(Q)^DRbar
  3  3     8.91343314E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.21193855E-05   # Y_d(Q)^DRbar
  2  2     1.72450643E-03   # Y_s(Q)^DRbar
  3  3     7.62749743E-02   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.74044671E-05   # Y_e(Q)^DRbar
  2  2     3.59867873E-03   # Y_mu(Q)^DRbar
  3  3     6.05215646E-02   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3     3.03370000E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     0.00000000E+00   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     0.00000000E+00   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    3.35000000E+02  # M_1
   2    1.95410000E+03  # M_2
   3    1.94560000E+03  # M_3
  21   -7.54590187E+06  # M^2_(H,d)
  22   -8.40774300E+06  # M^2_(H,u)
  31    3.26390000E+03  # M_(L,11)
  32    3.26390000E+03  # M_(L,22)
  33    3.26390000E+03  # M_(L,33)
  34    3.26390000E+03  # M_(E,11)
  35    3.26390000E+03  # M_(E,22)
  36    3.26390000E+03  # M_(E,33)
  41    3.21510000E+03  # M_(Q,11)
  42    3.21510000E+03  # M_(Q,22)
  43    1.57800000E+03  # M_(Q,33)
  44    3.21510000E+03  # M_(U,11)
  45    3.21510000E+03  # M_(U,22)
  46    6.09200000E+02  # M_(U,33)
  47    3.21510000E+03  # M_(D,11)
  48    3.21510000E+03  # M_(D,22)
  49    1.57800000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.03242013E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.21938365E+02  # h0
        35     9.42729415E+02  # H0
        36     9.48800000E+02  # A0
        37     9.51665722E+02  # H+
   1000001     3.22942065E+03  # ~d_L
   2000001     3.22767836E+03  # ~d_R
   1000002     3.22870127E+03  # ~u_L
   2000002     3.22814815E+03  # ~u_R
   1000003     3.22943026E+03  # ~s_L
   2000003     3.22766902E+03  # ~s_R
   1000004     3.22870574E+03  # ~c_L
   2000004     3.22814441E+03  # ~c_R
   1000005     1.51564867E+03  # ~b_1
   2000005     1.60674667E+03  # ~b_2
   1000006     3.75620823E+02  # ~t_1
   2000006     1.54180187E+03  # ~t_2
   1000011     3.26610288E+03  # ~e_L-
   2000011     3.26697967E+03  # ~e_R-
   1000012     3.26486183E+03  # ~nu_eL
   1000013     3.26602335E+03  # ~mu_L-
   2000013     3.26705948E+03  # ~mu_R-
   1000014     3.26486193E+03  # ~nu_muL
   1000015     3.26191828E+03  # ~tau_1-
   2000015     3.27124344E+03  # ~tau_2-
   1000016     3.26489071E+03  # ~nu_tauL
   1000021     2.08313434E+03  # ~g
   1000022     3.38277131E+02  # ~chi_10
   1000023     1.99653593E+03  # ~chi_20
   1000025    -2.90469810E+03  # ~chi_30
   1000035     2.90867080E+03  # ~chi_40
   1000024     1.99662369E+03  # ~chi_1+
   1000037     2.90898576E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -1.73456318E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    2.92580000E+03  # mu
   2    5.80000000E+00  # tan[beta](Q)
   3    2.44354501E+02  # v(Q)
   4    1.16939545E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1    -1.79533054E-01   # Re[R_st(1,1)]
   1  2     9.83751942E-01   # Re[R_st(1,2)]
   2  1    -9.83751942E-01   # Re[R_st(2,1)]
   2  2    -1.79533054E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.90044387E-01   # Re[R_sb(1,1)]
   1  2     1.40755502E-01   # Re[R_sb(1,2)]
   2  1    -1.40755502E-01   # Re[R_sb(2,1)]
   2  2     9.90044387E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     7.42909855E-01   # Re[R_sta(1,1)]
   1  2     6.69391475E-01   # Re[R_sta(1,2)]
   2  1    -6.69391475E-01   # Re[R_sta(2,1)]
   2  2     7.42909855E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99861162E-01   # Re[N(1,1)]
   1  2    -2.25293406E-04   # Re[N(1,2)]
   1  3     1.61413902E-02   # Re[N(1,3)]
   1  4    -4.13056359E-03   # Re[N(1,4)]
   2  1     1.25001739E-03   # Re[N(2,1)]
   2  2     9.97773111E-01   # Re[N(2,2)]
   2  3    -5.32216782E-02   # Re[N(2,3)]
   2  4     4.01834435E-02   # Re[N(2,4)]
   3  1    -8.49103532E-03   # Re[N(3,1)]
   3  2     9.24839770E-03   # Re[N(3,2)]
   3  3     7.07034894E-01   # Re[N(3,3)]
   3  4     7.07067202E-01   # Re[N(3,4)]
   4  1    -1.42827366E-02   # Re[N(4,1)]
   4  2     6.60547886E-02   # Re[N(4,2)]
   4  3     7.04988345E-01   # Re[N(4,3)]
   4  4    -7.05991644E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.96893025E-01   # Re[U(1,1)]
   1  2     7.87673569E-02   # Re[U(1,2)]
   2  1     7.87673569E-02   # Re[U(2,1)]
   2  2     9.96893025E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.98189979E-01   # Re[V(1,1)]
   1  2     6.01395560E-02   # Re[V(1,2)]
   2  1     6.01395560E-02   # Re[V(2,1)]
   2  2     9.98189979E-01   # Re[V(2,2)]
DECAY   1000011     1.94000298E+01   # ~e^-_L
#    BR                NDA      ID1      ID2
     2.13473911E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     2.62600804E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     5.23472266E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
     3.56060335E-04    2    -1000037        12   # BR(~e^-_L -> chi^-_2 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000011     1.65834484E+01   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99985081E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000013     1.91886881E+01   # ~mu^-_L
#    BR                NDA      ID1      ID2
     2.64958090E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     2.45415919E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     4.89221448E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     3.19015041E-04    2    -1000037        14   # BR(~mu^-_L -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.67948904E+01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.31259951E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     2.29397667E-02    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     4.57205536E-02    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     1.82175302E+01   # ~tau^-_1
#    BR                NDA      ID1      ID2
     5.33445877E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.55676628E-01    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     2.41307625E-04    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.20067965E-04    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     3.10514461E-01    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
DECAY   2000015     1.77976680E+01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     6.18453978E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     1.26898811E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     3.43977465E-04    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     5.73338712E-04    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     2.52750892E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     9.79002690E-04    2    -1000037        16   # BR(~tau^-_2 -> chi^-_2 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     1.93955827E+01   # ~nu_e
#    BR                NDA      ID1      ID2
     2.13774339E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     2.61575047E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     1.57299843E-04    2     1000035        12   # BR(~nu_e -> chi^0_4 nu_e)
     5.24281198E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
     2.06323835E-04    2     1000037        11   # BR(~nu_e -> chi^+_2 e^-)
DECAY   1000014     1.93956216E+01   # ~nu_mu
#    BR                NDA      ID1      ID2
     2.13773917E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     2.61574551E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     1.57299609E-04    2     1000035        14   # BR(~nu_mu -> chi^0_4 nu_mu)
     5.24280275E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     2.08165504E-04    2     1000037        13   # BR(~nu_mu -> chi^+_2 mu^-)
DECAY   1000016     1.94065853E+01   # ~nu_tau
#    BR                NDA      ID1      ID2
     2.13655111E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     2.61434586E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     1.57233459E-04    2     1000035        16   # BR(~nu_tau -> chi^0_4 nu_tau)
     5.24020280E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     7.26999347E-04    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000001     6.44993815E+01   # ~d_R
#    BR                NDA      ID1      ID2
     2.82085098E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.71788321E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     7.79351253E+01   # ~d_L
#    BR                NDA      ID1      ID2
     5.85378087E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     6.27635393E-02    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.25351799E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     8.05929372E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     6.45720775E+01   # ~s_R
#    BR                NDA      ID1      ID2
     2.80627703E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     4.14681787E-04    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     8.28387977E-04    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
     9.70693698E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     7.78624913E+01   # ~s_L
#    BR                NDA      ID1      ID2
     5.95377039E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     6.24789970E-02    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.24783173E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     8.06681900E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     1.06845316E+02   # ~b_1
#    BR                NDA      ID1      ID2
     1.98042301E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     5.20870123E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
     4.77149454E-01    2     1000006       -37   # BR(~b_1 -> ~t_1 H^-)
DECAY   2000005     3.31902364E+00   # ~b_2
#    BR                NDA      ID1      ID2
     2.50293086E-01    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     4.13115383E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
     3.36591531E-01    2     1000006       -37   # BR(~b_2 -> ~t_1 H^-)
DECAY   2000002     6.99938978E+01   # ~u_R
#    BR                NDA      ID1      ID2
     1.03993347E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     8.96005534E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     7.78951377E+01   # ~u_L
#    BR                NDA      ID1      ID2
     5.82751252E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     6.28077458E-02    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.25646939E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     8.05651535E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     7.00636142E+01   # ~c_R
#    BR                NDA      ID1      ID2
     1.03150440E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     5.34411607E-04    2     1000023         4   # BR(~c_R -> chi^0_2 c)
     1.06942405E-03    2     1000024         3   # BR(~c_R -> chi^+_1 s)
     8.95244427E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     7.78254609E+01   # ~c_L
#    BR                NDA      ID1      ID2
     6.49839790E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     6.23828660E-02    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.24796710E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     8.06254618E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     4.00337760E-10   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     1.09314223E+02   # ~t_2
#    BR                NDA      ID1      ID2
     2.87634591E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     2.73174363E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     2.57425349E-01    2     1000006        36   # BR(~t_2 -> ~t_1 A^0)
     2.36884724E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
     2.29639219E-01    2     1000006        35   # BR(~t_2 -> ~t_1 H^0)
DECAY   1000024     8.41774351E+00   # chi^+_1
#    BR                NDA      ID1      ID2
     2.43105314E-02    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     4.60017033E-01    2     2000006        -5   # BR(chi^+_1 -> ~t_2 b_bar)
     5.07630745E-01    2    -1000005         6   # BR(chi^+_1 -> ~b^*_1 t)
     6.05148128E-03    2    -2000005         6   # BR(chi^+_1 -> ~b^*_2 t)
     7.61287074E-04    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
     1.13094364E-03    2          37   1000022   # BR(chi^+_1 -> H^+ chi^0_1)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000037     1.23410851E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     5.29316593E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     7.50008175E-03    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     2.85052452E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     6.40457616E-03    2    -2000005         6   # BR(chi^+_2 -> ~b^*_2 t)
     1.62688416E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     4.85259339E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     1.09462818E-02    2          37   1000022   # BR(chi^+_2 -> H^+ chi^0_1)
     4.51764380E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     4.67852749E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.38760702E-03    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     1.36120903E-03    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     1.23826221E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     8.40230505E+00   # chi^0_2
#    BR                NDA      ID1      ID2
     1.43805636E-02    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     1.43805636E-02    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     2.31930897E-01    2     2000006        -6   # BR(chi^0_2 -> ~t_2 t_bar)
     2.31930897E-01    2    -2000006         6   # BR(chi^0_2 -> ~t^*_2 t)
     2.49619245E-01    2     1000005        -5   # BR(chi^0_2 -> ~b_1 b_bar)
     2.49619245E-01    2    -1000005         5   # BR(chi^0_2 -> ~b^*_1 b)
     3.05002926E-03    2     2000005        -5   # BR(chi^0_2 -> ~b_2 b_bar)
     3.05002926E-03    2    -2000005         5   # BR(chi^0_2 -> ~b^*_2 b)
     1.32798276E-04    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     4.06039304E-04    2     1000022        36   # BR(chi^0_2 -> chi^0_1 A^0)
     7.50604240E-04    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
     6.69380693E-04    2     1000022        35   # BR(chi^0_2 -> chi^0_1 H^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000025     1.22783947E+02   # chi^0_3
#    BR                NDA      ID1      ID2
     2.81395146E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     2.81395146E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     1.35930583E-01    2     2000006        -6   # BR(chi^0_3 -> ~t_2 t_bar)
     1.35930583E-01    2    -2000006         6   # BR(chi^0_3 -> ~t^*_2 t)
     1.06014020E-03    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     1.06014020E-03    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     1.02343169E-03    2     2000005        -5   # BR(chi^0_3 -> ~b_2 b_bar)
     1.02343169E-03    2    -2000005         5   # BR(chi^0_3 -> ~b^*_2 b)
     4.39613927E-02    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     4.39613927E-02    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     1.24998167E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     4.00301648E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     5.09102663E-03    2     1000022        36   # BR(chi^0_3 -> chi^0_1 A^0)
     3.88663900E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     7.78184395E-04    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
     6.08999039E-03    2     1000022        35   # BR(chi^0_3 -> chi^0_1 H^0)
#    BR                NDA      ID1      ID2       ID3
     1.55106821E-03    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     8.62769261E-04    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     1.22142710E-03    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     1.22142710E-03    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     1.20343051E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     2.68938971E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     2.68938971E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     1.51072392E-01    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     1.51072392E-01    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     1.26822328E-03    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     1.26822328E-03    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     1.18723253E-03    2     2000005        -5   # BR(chi^0_4 -> ~b_2 b_bar)
     1.18723253E-03    2    -2000005         5   # BR(chi^0_4 -> ~b^*_2 b)
     3.69154343E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     3.69154343E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     4.47167367E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     8.16541079E-04    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     6.12328571E-03    2     1000022        36   # BR(chi^0_4 -> chi^0_1 A^0)
     1.28003263E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     4.72127564E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
     5.11640818E-03    2     1000022        35   # BR(chi^0_4 -> chi^0_1 H^0)
#    BR                NDA      ID1      ID2       ID3
     1.45031971E-03    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     6.39822860E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     1.28274541E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     1.28274541E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     6.90691895E+01   # ~g
#    BR                NDA      ID1      ID2
     3.19989701E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     3.19989701E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     5.52632384E-02    2     2000006        -6   # BR(~g -> ~t_2 t_bar)
     5.52632384E-02    2    -2000006         6   # BR(~g -> ~t^*_2 t)
     7.12759277E-02    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     7.12759277E-02    2    -1000005         5   # BR(~g -> ~b^*_1 b)
     5.30787466E-02    2     2000005        -5   # BR(~g -> ~b_2 b_bar)
     5.30787466E-02    2    -2000005         5   # BR(~g -> ~b^*_2 b)
#    BR                NDA      ID1      ID2       ID3
     7.48014583E-04    3     1000022         6        -6   # BR(~g -> chi^0_1 t t_bar)
DECAY        25     2.87903636E-03   # h^0
#    BR                NDA      ID1      ID2
     3.38678699E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     9.56401747E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     3.02978144E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.78012896E-01    2           5        -5   # BR(h^0 -> b b_bar)
     4.22041182E-02    2           4        -4   # BR(h^0 -> c c_bar)
     9.85223686E-02    2          21        21   # BR(h^0 -> g g)
     3.45727877E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     8.24551239E-03    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     8.24551239E-03    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     8.24551239E-03    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     2.88592934E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     2.88592934E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     8.24551239E-03    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     8.24551239E-03    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     8.24551239E-03    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     2.88592934E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     2.88592934E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     5.81361505E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     5.81361505E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     4.98309861E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     3.32206574E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     2.32544602E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     2.32544602E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     2.32544602E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     2.32544602E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     2.32544602E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     2.20494712E+01   # H^0
#    BR                NDA      ID1      ID2
     3.02271602E-03    2          15       -15   # BR(H^0 -> tau^- tau^+)
     1.63568191E-02    2           5        -5   # BR(H^0 -> b b_bar)
     5.95712378E-02    2           6        -6   # BR(H^0 -> t t_bar)
     9.18980144E-01    2     1000006  -1000006   # BR(H^0 -> ~t_1 ~t^*_1)
     4.80920196E-04    2          23        23   # BR(H^0 -> Z Z)
     9.21231509E-04    2          24       -24   # BR(H^0 -> W^+ W^-)
     4.69715677E-04    2          25        25   # BR(H^0 -> h^0 h^0)
     1.67076572E-04    2          21        21   # BR(H^0 -> g g)
DECAY        36     1.86325676E+00   # A^0
#    BR                NDA      ID1      ID2
     1.27407902E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     3.60352440E-02    2          15       -15   # BR(A^0 -> tau^- tau^+)
     1.95009481E-01    2           5        -5   # BR(A^0 -> b b_bar)
     7.67225212E-01    2           6        -6   # BR(A^0 -> t t_bar)
     2.55588715E-04    2     1000022   1000022   # BR(A^0 -> chi^0_1 chi^0_1)
     9.96527880E-04    2          25        23   # BR(A^0 -> h^0 Z)
     2.35507711E-04    2          21        21   # BR(A^0 -> g g)
DECAY        37     1.58523805E+00   # H^+
#    BR                NDA      ID1      ID2
     1.50204979E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     4.24830268E-02    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     1.17508685E-04    2          -3         4   # BR(H^+ -> s_bar c)
     9.56091433E-01    2          -5         6   # BR(H^+ -> b_bar t)
     1.15758842E-03    2          25        24   # BR(H^+ -> h^0 W^+)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.03180636E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    3.36081936E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    3.36400000E+01        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99054508E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    3.06720083E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    2.97265161E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.03180636E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    3.36081936E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    3.36400000E+01        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99992596E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    7.40440305E-06        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99992596E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    7.40440305E-06        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15807800E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    9.44091705E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.85395025E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    7.40440305E-06        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99992596E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.24766806E-04   # BR(b -> s gamma)
    2    1.59515099E-06   # BR(b -> s mu+ mu-)
    3    3.54766605E-05   # BR(b -> s nu nu)
    4    5.18876715E-11   # BR(Bd -> mu+ mu-)
    5    3.38960483E-09   # BR(Bs -> mu+ mu-)
    6    1.10443130E-04   # BR(B_u -> tau nu)
    7    9.98676170E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73468693E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.02962805E+01   # |Delta(M_Bs)| [ps^-1] 
   20    8.94499472E-16   # Delta(g-2)_electron/2
   21    1.01945648E-11   # Delta(g-2)_muon/2
   22   -4.70599712E-09   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    4.40837558E-04   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87741858E-01   # C7
     0305 4422   00   1    -2.67330225E-02   # C7
     0305 4322   00   1    -6.15355155E-04   # C7'
     0305 6421   00   0    -9.48844135E-02   # C8
     0305 6421   00   1    -3.11958832E-02   # C8
     0305 6321   00   1    -7.15739134E-04   # C8'
 03051111 4133   00   0    -4.76154920E-01   # C9 e+e-
 03051111 4133   00   1     1.39287161E-04   # C9 e+e-
 03051111 4233   00   1    -1.26728804E-06   # C9' e+e-
 03051111 4137   00   0    -3.93086627E+00   # C10 e+e-
 03051111 4137   00   1    -5.60780326E-03   # C10 e+e-
 03051111 4237   00   1     3.72012735E-05   # C10' e+e-
 03051313 4133   00   0    -4.76154920E-01   # C9 mu+mu-
 03051313 4133   00   1     1.39286060E-04   # C9 mu+mu-
 03051313 4233   00   1    -1.26729635E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.93086627E+00   # C10 mu+mu-
 03051313 4137   00   1    -5.60780436E-03   # C10 mu+mu-
 03051313 4237   00   1     3.72012818E-05   # C10' mu+mu-
 03051212 4237   00   0     1.48992531E+00   # C11 nu_1 nu_1
 03051212 4237   00   1     1.35172330E-03   # C11 nu_1 nu_1
 03051212 4137   00   1    -8.97714464E-06   # C11' nu_1 nu_1
 03051414 4237   00   0     1.48992531E+00   # C11 nu_2 nu_2
 03051414 4237   00   1     1.35172330E-03   # C11 nu_2 nu_2
 03051414 4137   00   1    -8.97714264E-06   # C11' nu_2 nu_2
 03051616 4237   00   0     1.48992531E+00   # C11 nu_3 nu_3
 03051616 4237   00   1     1.35172330E-03   # C11 nu_3 nu_3
 03051616 4137   00   1    -8.97657947E-06   # C11' nu_3 nu_3
