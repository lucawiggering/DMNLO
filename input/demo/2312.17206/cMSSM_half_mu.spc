# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v4.0.5  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101;
# W. Porod, F.~Staub, Comput. Phys. Commun. 183 (2012) 2458
#                     arXiv:1104.1573 [hep-ph] 
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 23.10.2023,  14:28
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v4.0.5      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
Block MODSEL  # Model selection
    1    1    # mSUGRA model
Block MINPAR  # Input parameters
    1    3.00000000E+03  # m0      
    2    1.40000000E+03  # m12     
    3    2.00000000E+01  # tanb at m_Z   
    4    1.00000000E+00  # cos(phase_mu)
    5    1.20000000E+04  # A0
#
Block gauge Q=  1.22583636E+16  # (GUT scale)
   1    6.98162113E-01  # g'(M_GUT)^DRbar
   2    6.98162113E-01  # g(M_GUT)^DRbar
   3    6.96325898E-01  # g3(M_GUT)^DRbar
Block Yu Q=  1.22583636E+16  # (GUT scale)
  1  1     3.44688336E-06   # Y_u(M_GUT)^DRbar
  2  2     1.75102646E-03   # Y_c(M_GUT)^DRbar
  3  3     5.56477321E-01   # Y_t(M_GUT)^DRbar
Block Yd Q=  1.22583636E+16  # (GUT scale)
  1  1     9.66828437E-05   # Y_d(M_GUT)^DRbar
  2  2     1.83697538E-03   # Y_s(M_GUT)^DRbar
  3  3     9.67447274E-02   # Y_b(M_GUT)^DRbar
Block Ye Q=  1.22583636E+16  # (GUT scale)
  1  1     4.21350519E-05   # Y_e(M_GUT)^DRbar
  2  2     8.71267990E-03   # Y_mu(M_GUT)^DRbar
  3  3     1.48783905E-01   # Y_tau(M_GUT)^DRbar
Block SMINPUTS  # SM parameters
         1     1.27931861E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637900E-05  # G_mu [GeV^-2]
         3     1.18400000E-01  # alpha_s(MZ)^MSbar
         4     9.11876000E+01  # m_Z(pole)
         5     4.18000000E+00  # m_b(m_b), MSbar
         6     1.73100000E+02  # m_t(pole)
         7     1.77682000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998930E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658372E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     2.50000000E-03  # m_u(2 GeV), MSbar
        23     9.50000000E-02  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  1.36820083E+03  # (SUSY scale)
   1    3.61986117E-01  # g'(Q)^DRbar
   2    6.33800572E-01  # g(Q)^DRbar
   3    1.02128805E+00  # g3(Q)^DRbar
Block Yu Q=  1.36820083E+03  # (SUSY scale)
  1  1     7.03110679E-06   # Y_u(Q)^DRbar
  2  2     3.57180250E-03   # Y_c(Q)^DRbar
  3  3     8.56555990E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.36820083E+03  # (SUSY scale)
  1  1     2.51029069E-04   # Y_d(Q)^DRbar
  2  2     4.76952438E-03   # Y_s(Q)^DRbar
  3  3     2.26150023E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.36820083E+03  # (SUSY scale)
  1  1     5.89151066E-05   # Y_e(Q)^DRbar
  2  2     1.21816876E-02   # Y_mu(Q)^DRbar
  3  3     2.04307477E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.36820083E+03  # (SUSY scale)
  1  1     5.07357735E+03   # A_u(Q)^DRbar
  2  2     5.07348494E+03   # A_c(Q)^DRbar
  3  3     1.70956201E+03   # A_t(Q)^DRbar
Block Ad Q=  1.36820083E+03  # (SUSY scale)
  1  1     8.08275298E+03   # A_d(Q)^DRbar
  2  2     8.08262886E+03   # A_s(Q)^DRbar
  3  3     6.75931440E+03   # A_b(Q)^DRbar
Block Ae Q=  1.36820083E+03  # (SUSY scale)
  1  1     1.08162088E+04   # A_e(Q)^DRbar
  2  2     1.08147800E+04   # A_mu(Q)^DRbar
  3  3     1.04144932E+04   # A_tau(Q)^DRbar
Block MSOFT Q=  1.36820083E+03  # soft SUSY breaking masses at Q
   1    6.02248648E+02  # M_1
   2    1.08801335E+03  # M_2
   3    2.91587486E+03  # M_3
  21    5.17278006E+06  # M^2_(H,d)
  22   -1.26052812E+07  # M^2_(H,u)
  31    3.12173544E+03  # M_(L,11)
  32    3.12070090E+03  # M_(L,22)
  33    2.82571844E+03  # M_(L,33)
  34    3.03676161E+03  # M_(E,11)
  35    3.03462589E+03  # M_(E,22)
  36    2.38421605E+03  # M_(E,33)
  41    3.94762628E+03  # M_(Q,11)
  42    3.94754435E+03  # M_(Q,22)
  43    2.73112497E+03  # M_(Q,33)
  44    3.88291676E+03  # M_(U,11)
  45    3.88286705E+03  # M_(U,22)
  46    6.74093262E+02  # M_(U,33)
  47    3.87406357E+03  # M_(D,11)
  48    3.87394540E+03  # M_(D,22)
  49    3.62807508E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73100000E+02  # m_t(pole)
        23     9.11876000E+01  # m_Z(pole)
        24     8.03588098E+01  # W+
        15     1.77682000E+00  # m_tau(pole)
        25     1.18904933E+02  # h0
        35     4.22791686E+03  # H0
        36     4.22771963E+03  # A0
        37     4.22887288E+03  # H+
   1000001     3.96534016E+03  # ~d_L
   2000001     3.88668328E+03  # ~d_R
   1000002     3.96468285E+03  # ~u_L
   2000002     3.89699812E+03  # ~u_R
   1000003     3.96526332E+03  # ~s_L
   2000003     3.88657104E+03  # ~s_R
   1000004     3.96460676E+03  # ~c_L
   2000004     3.89694771E+03  # ~c_R
   1000005     2.76590845E+03  # ~b_1
   2000005     3.65124598E+03  # ~b_2
   1000006     6.48338350E+02  # ~t_1
   2000006     2.77099037E+03  # ~t_2
   1000011     3.12778340E+03  # ~e_L-
   2000011     3.03835480E+03  # ~e_R-
   1000012     3.12645274E+03  # ~nu_eL
   1000013     3.12679429E+03  # ~mu_L-
   2000013     3.03627908E+03  # ~mu_R-
   1000014     3.12545181E+03  # ~nu_muL
   1000015     2.40778240E+03  # ~tau_1-
   2000015     2.84103580E+03  # ~tau_2-
   1000016     2.83885250E+03  # ~nu_tauL
   1000021     3.05723004E+03  # ~g
   1000022     6.06294841E+02  # ~chi_10
   1000023     1.14339310E+03  # ~chi_20
   1000025    -3.55600029E+03  # ~chi_30
   1000035     3.55690583E+03  # ~chi_40
   1000024     1.14355557E+03  # ~chi_1+
   1000037     3.55751419E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -5.19976029E-02   # alpha
Block Hmix Q=  1.36820083E+03  # Higgs mixing parameters
   1    3.56189106E+03  # mu
   2    1.92609098E+01  # tan[beta](Q)
   3    2.44187426E+02  # v(Q)
   4    1.78662220E+07  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1    -3.34800646E-02   # Re[R_st(1,1)]
   1  2     9.99439385E-01   # Re[R_st(1,2)]
   2  1    -9.99439385E-01   # Re[R_st(2,1)]
   2  2    -3.34800646E-02   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.99777014E-01   # Re[R_sb(1,1)]
   1  2     2.11168630E-02   # Re[R_sb(1,2)]
   2  1    -2.11168630E-02   # Re[R_sb(2,1)]
   2  2     9.99777014E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     4.71836105E-02   # Re[R_sta(1,1)]
   1  2     9.98886233E-01   # Re[R_sta(1,2)]
   2  1    -9.98886233E-01   # Re[R_sta(2,1)]
   2  2     4.71836105E-02   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99902574E-01   # Re[N(1,1)]
   1  2    -4.41104089E-04   # Re[N(1,2)]
   1  3     1.36438991E-02   # Re[N(1,3)]
   1  4    -2.91410237E-03   # Re[N(1,4)]
   2  1    -7.98700574E-04   # Re[N(2,1)]
   2  2    -9.99665280E-01   # Re[N(2,2)]
   2  3     2.43631184E-02   # Re[N(2,3)]
   2  4    -8.66763901E-03   # Re[N(2,4)]
   3  1     7.58111688E-03   # Re[N(3,1)]
   3  2    -1.11048556E-02   # Re[N(3,2)]
   3  3    -7.06973435E-01   # Re[N(3,3)]
   3  4    -7.07112276E-01   # Re[N(3,4)]
   4  1    -1.16932157E-02   # Re[N(4,1)]
   4  2     2.33626847E-02   # Re[N(4,2)]
   4  3     7.06688648E-01   # Re[N(4,3)]
   4  4    -7.07042155E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.99405857E-01   # Re[U(1,1)]
   1  2     3.44663938E-02   # Re[U(1,2)]
   2  1     3.44663938E-02   # Re[U(2,1)]
   2  2     9.99405857E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.99924543E-01   # Re[V(1,1)]
   1  2     1.22844444E-02   # Re[V(1,2)]
   2  1     1.22844444E-02   # Re[V(2,1)]
   2  2     9.99924543E-01   # Re[V(2,2)]
###########################################################
# start writing paramaters at Q=  9.11876000E+01
###########################################################
Block gauge Q=  9.11876000E+01  # (SUSY scale)
   1    3.53205237E-01  # g'(Q)^DRbar
   2    6.28620341E-01  # g(Q)^DRbar
   3    1.07873187E+00  # g3(Q)^DRbar
Block Yu Q=  9.11876000E+01  # (SUSY scale)
  1  1     7.64471982E-06   # Y_u(Q)^DRbar
  2  2     3.88351349E-03   # Y_c(Q)^DRbar
  3  3     8.94995210E-01   # Y_t(Q)^DRbar
Block Yd Q=  9.11876000E+01  # (SUSY scale)
  1  1     2.82712360E-04   # Y_d(Q)^DRbar
  2  2     5.37149502E-03   # Y_s(Q)^DRbar
  3  3     2.50710954E-01   # Y_b(Q)^DRbar
Block Ye Q=  9.11876000E+01  # (SUSY scale)
  1  1     6.02930857E-05   # Y_e(Q)^DRbar
  2  2     1.24665101E-02   # Y_mu(Q)^DRbar
  3  3     2.08626362E-01   # Y_tau(Q)^DRbar
Block EXTPAR  # non-universal input parameters
   0    9.11876000E+01  # Q
  25    1.99999921E+01  # tan(beta)
Block Au Q=  9.11876000E+01  # (SUSY scale)
  1  1     4.31607060E+03   # A_u(Q)^DRbar
  2  2     4.31596483E+03   # A_c(Q)^DRbar
  3  3     8.43713262E+02   # A_t(Q)^DRbar
Block Ad Q=  9.11876000E+01  # (SUSY scale)
  1  1     7.36562734E+03   # A_d(Q)^DRbar
  2  2     7.36548067E+03   # A_s(Q)^DRbar
  3  3     5.97334325E+03   # A_b(Q)^DRbar
Block Ae Q=  9.11876000E+01  # (SUSY scale)
  1  1     1.07124496E+04   # A_e(Q)^DRbar
  2  2     1.07108522E+04   # A_mu(Q)^DRbar
  3  3     1.02653001E+04   # A_tau(Q)^DRbar
Block MSOFT Q=  9.11876000E+01  # soft SUSY breaking masses at Q
   1    5.70197884E+02  # M_1
   2    1.06162398E+03  # M_2
   3    3.23873997E+03  # M_3
  21    4.62853829E+06  # M^2_(H,d)
  22   -1.24997326E+07  # M^2_(H,u)
  31    3.12834088E+03  # M_(L,11)
  32    3.12719239E+03  # M_(L,22)
  33    2.80064586E+03  # M_(L,33)
  34    3.03717545E+03  # M_(E,11)
  35    3.03479967E+03  # M_(E,22)
  36    2.30660054E+03  # M_(E,33)
  41    4.17696941E+03  # M_(Q,11)
  42    4.17687948E+03  # M_(Q,22)
  43    3.03693419E+03  # M_(Q,33)
  44    4.11370520E+03  # M_(U,11)
  45    4.11365348E+03  # M_(U,22)
  46    1.53773080E+03  # M_(U,33)
  47    4.10524552E+03  # M_(D,11)
  48    4.10511323E+03  # M_(D,22)
  49    3.83999874E+03  # M_(D,33)
Block Hmix Q=  9.11876000E+01  # Higgs mixing parameters
   1    3.48596219E+03  # mu
   2    1.99999921E+01  # tan[beta](Q)
   4    2.68823920E+06  # m^2_A(Q)
###########################################################
# stop writing paramaters at Q=  9.11876000E+01
###########################################################
###########################################################
# start writing paramaters at Q=  6.84100000E+02
###########################################################
Block gauge Q=  6.84100000E+02  # (SUSY scale)
   1    3.59677874E-01  # g'(Q)^DRbar
   2    6.32468328E-01  # g(Q)^DRbar
   3    1.03512644E+00  # g3(Q)^DRbar
Block Yu Q=  6.84100000E+02  # (SUSY scale)
  1  1     7.17867429E-06   # Y_u(Q)^DRbar
  2  2     3.64676580E-03   # Y_c(Q)^DRbar
  3  3     8.65963291E-01   # Y_t(Q)^DRbar
Block Yd Q=  6.84100000E+02  # (SUSY scale)
  1  1     2.58541904E-04   # Y_d(Q)^DRbar
  2  2     4.91226562E-03   # Y_s(Q)^DRbar
  3  3     2.32017630E-01   # Y_b(Q)^DRbar
Block Ye Q=  6.84100000E+02  # (SUSY scale)
  1  1     5.92720668E-05   # Y_e(Q)^DRbar
  2  2     1.22554710E-02   # Y_mu(Q)^DRbar
  3  3     2.05431388E-01   # Y_tau(Q)^DRbar
Block EXTPAR  # non-universal input parameters
   0    6.84100000E+02  # Q
  25    1.94416795E+01  # tan(beta)
Block Au Q=  6.84100000E+02  # (SUSY scale)
  1  1     4.88599038E+03   # A_u(Q)^DRbar
  2  2     4.88589465E+03   # A_c(Q)^DRbar
  3  3     1.48893491E+03   # A_t(Q)^DRbar
Block Ad Q=  6.84100000E+02  # (SUSY scale)
  1  1     7.91171103E+03   # A_d(Q)^DRbar
  2  2     7.91158142E+03   # A_s(Q)^DRbar
  3  3     6.56916021E+03   # A_b(Q)^DRbar
Block Ae Q=  6.84100000E+02  # (SUSY scale)
  1  1     1.07897414E+04   # A_e(Q)^DRbar
  2  2     1.07882699E+04   # A_mu(Q)^DRbar
  3  3     1.03765155E+04   # A_tau(Q)^DRbar
Block MSOFT Q=  6.84100000E+02  # soft SUSY breaking masses at Q
   1    5.93788682E+02  # M_1
   2    1.08130626E+03  # M_2
   3    2.99220999E+03  # M_3
  21    5.03787903E+06  # M^2_(H,d)
  22   -1.25650993E+07  # M^2_(H,u)
  31    3.12348374E+03  # M_(L,11)
  32    3.12242030E+03  # M_(L,22)
  33    2.81941169E+03  # M_(L,33)
  34    3.03688639E+03  # M_(E,11)
  35    3.03468985E+03  # M_(E,22)
  36    2.36468216E+03  # M_(E,33)
  41    4.00089321E+03  # M_(Q,11)
  42    4.00080920E+03  # M_(Q,22)
  43    2.80383009E+03  # M_(Q,33)
  44    3.93638542E+03  # M_(U,11)
  45    3.93633514E+03  # M_(U,22)
  46    9.46369114E+02  # M_(U,33)
  47    3.92761757E+03  # M_(D,11)
  48    3.92749577E+03  # M_(D,22)
  49    3.67642546E+03  # M_(D,33)
Block Hmix Q=  6.84100000E+02  # Higgs mixing parameters
   1    3.54388117E+03  # mu
   2    1.94416795E+01  # tan[beta](Q)
   4    1.20574479E+07  # m^2_A(Q)
###########################################################
# stop writing paramaters at Q=  6.84100000E+02
###########################################################
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   2000011     1.46017273E+01   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99999484E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000011     3.18906878E+01   # ~e^-_L
#    BR                NDA      ID1      ID2
     1.18204805E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     2.94229501E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     5.87565694E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.45926566E+01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99729347E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     1.80016127E-04    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
DECAY   1000013     3.18722428E+01   # ~mu^-_L
#    BR                NDA      ID1      ID2
     1.18274932E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     2.94206177E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     5.87518828E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     1.10534608E+01   # ~tau^-_1
#    BR                NDA      ID1      ID2
     9.94683372E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.77368546E-03    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     3.54294295E-03    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
DECAY   2000015     2.79043592E+01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     1.21375149E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     2.84877064E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     5.68813449E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     1.21484577E-02    2     1000015        23   # BR(~tau^-_2 -> ~tau^-_1 Z)
     1.27858801E-02    2     1000015        25   # BR(~tau^-_2 -> ~tau^-_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     3.18834866E+01   # ~nu_e
#    BR                NDA      ID1      ID2
     1.18539020E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     2.93557179E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     5.87903801E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
DECAY   1000014     3.18675749E+01   # ~nu_mu
#    BR                NDA      ID1      ID2
     1.18554302E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     2.93551648E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     5.87892919E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
DECAY   1000016     2.79054775E+01   # ~nu_tau
#    BR                NDA      ID1      ID2
     1.20940851E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     2.84859761E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     5.70535014E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     2.36643745E-02    2     1000015        24   # BR(~nu_tau -> ~tau^-_1 W^+)
DECAY   2000001     3.34073157E+01   # ~d_R
#    BR                NDA      ID1      ID2
     6.41441328E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.35855418E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     7.65576641E+01   # ~d_L
#    BR                NDA      ID1      ID2
     7.18641116E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     1.73810748E-01    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     3.47528902E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     4.71448793E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     3.34013344E+01   # ~s_R
#    BR                NDA      ID1      ID2
     6.41530983E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     9.35822317E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     7.65511217E+01   # ~s_L
#    BR                NDA      ID1      ID2
     7.18709583E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     1.73817824E-01    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     3.47542922E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     4.71425549E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     3.52800557E+01   # ~b_1
#    BR                NDA      ID1      ID2
     1.04079938E-02    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     2.15175891E-01    2     1000023         5   # BR(~b_1 -> chi^0_2 b)
     4.24492742E-01    2    -1000024         6   # BR(~b_1 -> chi^-_1 t)
     3.49923373E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     2.22105999E+01   # ~b_2
#    BR                NDA      ID1      ID2
     8.99304056E-02    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     2.21259678E-04    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     2.18571455E-04    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     8.12494760E-01    2     1000021         5   # BR(~b_2 -> ~g b)
     6.14748235E-04    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
     4.54734607E-02    2     2000006       -24   # BR(~b_2 -> ~t_2 W^-)
     2.37010656E-02    2     1000005        23   # BR(~b_2 -> ~b_1 Z)
     2.72216576E-02    2     1000005        25   # BR(~b_2 -> ~b_1 h^0)
DECAY   2000002     4.04841676E+01   # ~u_R
#    BR                NDA      ID1      ID2
     2.12343450E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     7.87655222E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     7.65367224E+01   # ~u_L
#    BR                NDA      ID1      ID2
     7.12076302E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     1.73924829E-01    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     3.47906319E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     4.71041455E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     4.04813506E+01   # ~c_R
#    BR                NDA      ID1      ID2
     2.12349224E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     7.87620679E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     7.65305982E+01   # ~c_L
#    BR                NDA      ID1      ID2
     7.12437360E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     1.73929200E-01    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     3.47915190E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     4.71023225E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     9.77651415E-07   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     3.52554277E+01   # ~t_2
#    BR                NDA      ID1      ID2
     1.05177380E-02    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     2.13305091E-01    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     4.32365149E-01    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     1.81158265E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     1.62653756E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     2.08320925E-03   # chi^+_1
#    BR                NDA      ID1      ID2
     8.64796099E-01    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     1.30994298E-01    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     1.98422377E-03    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
     5.89240598E-04    3     1000022       -11        12   # BR(chi^+_1 -> chi^0_1 e^+ nu_e)
     5.90066071E-04    3     1000022       -13        14   # BR(chi^+_1 -> chi^0_1 mu^+ nu_mu)
     9.03716537E-04    3     1000022       -15        16   # BR(chi^+_1 -> chi^0_1 tau^+ nu_tau)
DECAY   1000037     1.10557828E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     3.87089733E-03    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     1.76269722E-03    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     6.58079828E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     7.65206641E-03    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     1.07071835E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     2.23460151E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     6.08164051E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     6.24972790E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     6.71544539E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     4.17656341E-03    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     1.77404706E-03    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     2.72512887E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     6.38935754E-03   # chi^0_2
#    BR                NDA      ID1      ID2
     3.28113025E-01    2     1000006        -4   # BR(chi^0_2 -> ~t_1 c_bar)
     3.28113025E-01    2    -1000006         4   # BR(chi^0_2 -> ~t^*_1 c)
     1.44943503E-01    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     1.44943503E-01    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     2.71162179E-03    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     5.02696330E-02    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.25245156E-04    3     1000022         6        -6   # BR(chi^0_2 -> chi^0_1 t t_bar)
     1.47505820E-04    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
     3.39097674E-04    3     1000022        12       -12   # BR(chi^0_2 -> chi^0_1 nu_e nu_bar_e)
DECAY   1000025     1.10058885E+02   # chi^0_3
#    BR                NDA      ID1      ID2
     1.95872398E-03    2     1000015       -15   # BR(chi^0_3 -> ~tau^-_1 tau^+)
     1.95872398E-03    2    -1000015        15   # BR(chi^0_3 -> ~tau^+_1 tau^-)
     8.81860854E-04    2     2000015       -15   # BR(chi^0_3 -> ~tau^-_2 tau^+)
     8.81860854E-04    2    -2000015        15   # BR(chi^0_3 -> ~tau^+_2 tau^-)
     3.32218142E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     3.32218142E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     5.24991836E-02    2     2000006        -6   # BR(chi^0_3 -> ~t_2 t_bar)
     5.24991836E-02    2    -2000006         6   # BR(chi^0_3 -> ~t^*_2 t)
     3.84391332E-03    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     3.84391332E-03    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     6.09986703E-02    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     6.09986703E-02    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     1.63183876E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     5.10961516E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     6.29615927E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     1.18779673E-02    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.91975295E-03    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     2.20305632E-03    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     2.73119275E-03    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     2.73119275E-03    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     1.09920391E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     1.95553917E-03    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     1.95553917E-03    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     8.88987337E-04    2     2000015       -15   # BR(chi^0_4 -> ~tau^-_2 tau^+)
     8.88987337E-04    2    -2000015        15   # BR(chi^0_4 -> ~tau^+_2 tau^-)
     3.28816547E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     3.28816547E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     5.40887389E-02    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     5.40887389E-02    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     3.85807566E-03    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     3.85807566E-03    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     6.09498392E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     6.09498392E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     6.86853277E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     1.15821588E-02    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     1.54594561E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     5.56016714E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.86422654E-03    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     1.99392289E-03    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     2.73162098E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     2.73162098E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     6.21139404E+01   # ~g
#    BR                NDA      ID1      ID2
     3.09631348E-04    2     1000006        -4   # BR(~g -> ~t_1 c_bar)
     3.09631348E-04    2    -1000006         4   # BR(~g -> ~t^*_1 c)
     4.69294449E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     4.69294449E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     1.26338642E-02    2     2000006        -6   # BR(~g -> ~t_2 t_bar)
     1.26338642E-02    2    -2000006         6   # BR(~g -> ~t^*_2 t)
     1.68099363E-02    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     1.68099363E-02    2    -1000005         5   # BR(~g -> ~b^*_1 b)
#    BR                NDA      ID1      ID2       ID3
     1.25183328E-03    3     1000022         6        -6   # BR(~g -> chi^0_1 t t_bar)
DECAY        25     2.33464761E-03   # h^0
#    BR                NDA      ID1      ID2
     4.23501935E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     1.18498834E-01    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.57578721E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.61013935E-01    2           5        -5   # BR(h^0 -> b b_bar)
     5.00546499E-02    2           4        -4   # BR(h^0 -> c c_bar)
     1.14039851E-01    2          21        21   # BR(h^0 -> g g)
     3.81575076E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     6.96718426E-03    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     6.96718426E-03    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     6.96718426E-03    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     2.43851449E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     2.43851449E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     6.96718426E-03    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     6.96718426E-03    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     6.96718426E-03    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     2.43851449E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     2.43851449E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     4.21723541E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     4.21723541E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     4.24233800E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     2.51025917E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     1.95800216E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     1.95800216E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     1.89524568E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     1.45595032E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     1.50615550E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     8.92995595E+01   # H^0
#    BR                NDA      ID1      ID2
     1.39395075E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     3.92103369E-02    2          15       -15   # BR(H^0 -> tau^- tau^+)
     1.66092950E-01    2           5        -5   # BR(H^0 -> b b_bar)
     6.43204974E-03    2           6        -6   # BR(H^0 -> t t_bar)
     2.92495016E-03    2     1000006  -1000006   # BR(H^0 -> ~t_1 ~t^*_1)
     3.90357643E-01    2     1000006  -2000006   # BR(H^0 -> ~t_1 ~t^*_2)
     3.90357643E-01    2    -1000006   2000006   # BR(H^0 -> ~t^*_1 ~t_2)
     1.04056045E-04    2     1000022   1000023   # BR(H^0 -> chi^0_1 chi^0_2)
     2.19790582E-03    2     1000022   1000025   # BR(H^0 -> chi^0_1 chi^0_3)
     1.04657176E-04    2     1000022   1000035   # BR(H^0 -> chi^0_1 chi^0_4)
     1.28182187E-04    2     1000023   1000023   # BR(H^0 -> chi^0_2 chi^0_2)
     2.56848685E-04    2     1000024  -1000024   # BR(H^0 -> chi^+_1 chi^-_1)
     5.52639396E-04    2          23        23   # BR(H^0 -> Z Z)
     1.04264483E-03    2          24       -24   # BR(H^0 -> W^+ W^-)
DECAY        36     8.92127031E+01   # A^0
#    BR                NDA      ID1      ID2
     1.39526101E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     3.92472208E-02    2          15       -15   # BR(A^0 -> tau^- tau^+)
     1.66248831E-01    2           5        -5   # BR(A^0 -> b b_bar)
     6.42595657E-03    2           6        -6   # BR(A^0 -> t t_bar)
     3.92418090E-01    2     1000006  -2000006   # BR(A^0 -> ~t_1 ~t^*_2)
     3.92418090E-01    2    -1000006   2000006   # BR(A^0 -> ~t^*_1 ~t_2)
     1.31209486E-04    2     1000022   1000023   # BR(A^0 -> chi^0_1 chi^0_2)
     1.06725388E-04    2     1000022   1000025   # BR(A^0 -> chi^0_1 chi^0_3)
     2.17508089E-03    2     1000022   1000035   # BR(A^0 -> chi^0_1 chi^0_4)
     1.95281706E-04    2     1000023   1000023   # BR(A^0 -> chi^0_2 chi^0_2)
     3.91385035E-04    2     1000024  -1000024   # BR(A^0 -> chi^+_1 chi^-_1)
DECAY        37     8.72322101E+01   # H^+
#    BR                NDA      ID1      ID2
     1.42732782E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     4.01492265E-02    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     1.52767106E-01    2          -5         6   # BR(H^+ -> b_bar t)
     8.04295620E-01    2    -1000005   1000006   # BR(H^+ -> ~b^*_1 ~t_1)
     2.40548835E-04    2     1000024   1000022   # BR(H^+ -> chi^+_1 chi^0_1)
     2.33462591E-03    2     1000037   1000022   # BR(H^+ -> chi^+_2 chi^0_1)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00484224E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    3.70977804E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    3.70982646E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99986948E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    2.70859634E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    2.69554387E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00484224E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    3.70977804E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    3.70982646E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999984E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    1.57627505E-08        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999984E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    1.57627505E-08        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.18624221E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    2.29350427E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    9.54325109E-02        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    1.57627505E-08        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999984E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.23272537E-04   # BR(b -> s gamma)
    2    1.58998104E-06   # BR(b -> s mu+ mu-)
    3    3.53262923E-05   # BR(b -> s nu nu)
    4    2.19301828E-15   # BR(Bd -> e+ e-)
    5    9.36832995E-11   # BR(Bd -> mu+ mu-)
    6    1.96126885E-08   # BR(Bd -> tau+ tau-)
    7    7.38816405E-14   # BR(Bs -> e+ e-)
    8    3.15622301E-09   # BR(Bs -> mu+ mu-)
    9    6.69497702E-07   # BR(Bs -> tau+ tau-)
   10    9.67224093E-05   # BR(B_u -> tau nu)
   11    9.99255652E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41414682E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93451573E+01   # |Delta(M_Bs)| [ps^-1] 
   14    0.00000000E+00   # neutron EDM according to the chiral quark model
   15    0.00000000E+00   # neutron EDM according to the relativistic quark-parton model
   16    2.13900762E-03   # epsilon_K
   17    2.93789973E-15   # Delta(M_K)
   18    2.48514893E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.30266726E-11   # BR(K^+ -> pi^+ nu nu)
   20    8.24355379E-16   # Delta(g-2)_electron/2
   21    3.52594439E-11   # Delta(g-2)_muon/2
   22    1.14400338E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    4.24276683E-05   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88923374E-01   # C7
     0305 4422   00   2    -1.92272110E-01   # C7
     0305 4322   00   2    -6.70015160E-05   # C7'
     0305 6421   00   0    -9.52609703E-02   # C8
     0305 6421   00   2    -9.91136062E-02   # C8
     0305 6321   00   2    -7.81046535E-05   # C8'
 03051111 4133   00   0     1.16716179E+00   # C9 e+e-
 03051111 4133   00   2     1.16725228E+00   # C9 e+e-
 03051111 4233   00   2    -5.32030851E-07   # C9' e+e-
 03051111 4137   00   0    -3.99512308E+00   # C10 e+e-
 03051111 4137   00   2    -3.99390051E+00   # C10 e+e-
 03051111 4237   00   2     1.51477139E-05   # C10' e+e-
 03051313 4133   00   0     1.16716179E+00   # C9 mu+mu-
 03051313 4133   00   2     1.16725199E+00   # C9 mu+mu-
 03051313 4233   00   2    -5.32036651E-07   # C9' mu+mu-
 03051313 4137   00   0    -3.99512308E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.99390080E+00   # C10 mu+mu-
 03051313 4237   00   2     1.51477162E-05   # C10' mu+mu-
 03051212 4137   00   0     1.50685489E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50656056E+00   # C11 nu_1 nu_1
 03051212 4237   00   2    -3.65310774E-06   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50685489E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50656056E+00   # C11 nu_2 nu_2
 03051414 4237   00   2    -3.65310667E-06   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50685489E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50656056E+00   # C11 nu_3 nu_3
 03051616 4237   00   2    -3.65280540E-06   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.86484367E-07   # C7
     0305 4422   00   2     9.99437973E-07   # C7
     0305 4322   00   2    -3.81048706E-07   # C7'
     0305 6421   00   0     3.31049081E-07   # C8
     0305 6421   00   2    -3.05551751E-08   # C8
     0305 6321   00   2    -2.89847363E-07   # C8'
 03051111 4133   00   2    -2.52879760E-08   # C9 e+e-
 03051111 4233   00   2     1.43491760E-08   # C9' e+e-
 03051111 4137   00   2     1.56654821E-06   # C10 e+e-
 03051111 4237   00   2    -4.03653678E-07   # C10' e+e-
 03051313 4133   00   2    -2.53003817E-08   # C9 mu+mu-
 03051313 4233   00   2     1.43491542E-08   # C9' mu+mu-
 03051313 4137   00   2     1.56656228E-06   # C10 mu+mu-
 03051313 4237   00   2    -4.03653718E-07   # C10' mu+mu-
 03051212 4137   00   2    -3.67583196E-07   # C11 nu_1 nu_1
 03051212 4237   00   2     9.73408898E-08   # C11' nu_1 nu_1
 03051414 4137   00   2    -3.67582848E-07   # C11 nu_2 nu_2
 03051414 4237   00   2     9.73408896E-08   # C11' nu_2 nu_2
 03051616 4137   00   2    -3.67483019E-07   # C11 nu_3 nu_3
 03051616 4237   00   2     9.73408427E-08   # C11' nu_3 nu_3
