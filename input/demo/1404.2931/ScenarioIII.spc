# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 09.06.2023,  11:33
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
    3    1.03603328E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.84800000E+03  # M_1
    2    1.80000000E+03  # M_2
    3    1.10230000E+03  # M_3
   23    1.10000000E+03  # mu
   25    1.00000000E+01  # tan(beta)
   26    1.95140000E+03  # m_A, pole mass
   31    1.98210000E+03  # M^2_L11
   32    1.98210000E+03  # M^2_L22
   33    1.98210000E+03  # M^2_L33
   34    1.98210000E+03  # M^2_E11
   35    1.98210000E+03  # M^2_E22
   36    1.98210000E+03  # M^2_E33
   41    3.98850000E+03  # M^2_Q11
   42    3.98850000E+03  # M^2_Q22
   43    2.30200000E+03  # M^2_Q33
   44    3.98850000E+03  # M^2_U11
   45    3.98850000E+03  # M^2_U22
   46    1.63660000E+03  # M^2_U33
   47    3.98850000E+03  # M^2_D11
   48    3.98850000E+03  # M^2_D22
   49    2.30200000E+03  # M^2_D33
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
   1    3.61312959E-01  # g'(Q)^DRbar
   2    6.33594463E-01  # g(Q)^DRbar
   3    1.04206005E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.61277098E-06   # Y_u(Q)^DRbar
  2  2     3.64607197E-03   # Y_c(Q)^DRbar
  3  3     8.57981806E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.43212386E-04   # Y_d(Q)^DRbar
  2  2     3.00745775E-03   # Y_s(Q)^DRbar
  3  3     1.43373483E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.99106981E-05   # Y_e(Q)^DRbar
  2  2     6.18455899E-03   # Y_mu(Q)^DRbar
  3  3     1.03999516E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3    -2.90833673E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     0.00000000E+00   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     0.00000000E+00   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    1.84800000E+03  # M_1
   2    1.80000000E+03  # M_2
   3    1.10230000E+03  # M_3
  21    2.52015371E+06  # M^2_(H,d)
  22   -1.38632624E+06  # M^2_(H,u)
  31    1.98210000E+03  # M_(L,11)
  32    1.98210000E+03  # M_(L,22)
  33    1.98210000E+03  # M_(L,33)
  34    1.98210000E+03  # M_(E,11)
  35    1.98210000E+03  # M_(E,22)
  36    1.98210000E+03  # M_(E,33)
  41    3.98850000E+03  # M_(Q,11)
  42    3.98850000E+03  # M_(Q,22)
  43    2.30200000E+03  # M_(Q,33)
  44    3.98850000E+03  # M_(U,11)
  45    3.98850000E+03  # M_(U,22)
  46    1.63660000E+03  # M_(U,33)
  47    3.98850000E+03  # M_(D,11)
  48    3.98850000E+03  # M_(D,22)
  49    2.30200000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02981454E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.26000168E+02  # h0
        35     1.95146776E+03  # H0
        36     1.95140000E+03  # A0
        37     1.95319665E+03  # H+
   1000001     4.02551727E+03  # ~d_L
   2000001     4.02185752E+03  # ~d_R
   1000002     4.02491302E+03  # ~u_L
   2000002     4.02155628E+03  # ~u_R
   1000003     4.02551930E+03  # ~s_L
   2000003     4.02185692E+03  # ~s_R
   1000004     4.02491389E+03  # ~c_L
   2000004     4.02155725E+03  # ~c_R
   1000005     2.31178977E+03  # ~b_1
   2000005     2.32644281E+03  # ~b_2
   1000006     1.66408447E+03  # ~t_1
   2000006     2.33491525E+03  # ~t_2
   1000011     1.98682392E+03  # ~e_L-
   2000011     1.98581137E+03  # ~e_R-
   1000012     1.98491914E+03  # ~nu_eL
   1000013     1.98690510E+03  # ~mu_L-
   2000013     1.98573015E+03  # ~mu_R-
   1000014     1.98491915E+03  # ~nu_muL
   1000015     1.98127550E+03  # ~tau_1-
   2000015     1.99134778E+03  # ~tau_2-
   1000016     1.98492153E+03  # ~nu_tauL
   1000021     1.30298582E+03  # ~g
   1000022     1.10670778E+03  # ~chi_10
   1000023    -1.11492905E+03  # ~chi_20
   1000025     1.85500149E+03  # ~chi_30
   1000035     1.86563106E+03  # ~chi_40
   1000024     1.10957660E+03  # ~chi_1+
   1000037     1.85632541E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -1.00251181E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    1.10000000E+03  # mu
   2    1.00000000E+01  # tan[beta](Q)
   3    2.45069263E+02  # v(Q)
   4    3.82386740E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     1.72507757E-01   # Re[R_st(1,1)]
   1  2     9.85008159E-01   # Re[R_st(1,2)]
   2  1    -9.85008159E-01   # Re[R_st(2,1)]
   2  2     1.72507757E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.13544831E-01   # Re[R_sb(1,1)]
   1  2     4.06738051E-01   # Re[R_sb(1,2)]
   2  1    -4.06738051E-01   # Re[R_sb(2,1)]
   2  2     9.13544831E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.72514178E-01   # Re[R_sta(1,1)]
   1  2     7.40084238E-01   # Re[R_sta(1,2)]
   2  1    -7.40084238E-01   # Re[R_sta(2,1)]
   2  2     6.72514178E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     4.60158606E-02   # Re[N(1,1)]
   1  2    -8.15954061E-02   # Re[N(1,2)]
   1  3     7.05779479E-01   # Re[N(1,3)]
   1  4    -7.02210835E-01   # Re[N(1,4)]
   2  1     9.59645265E-03   # Re[N(2,1)]
   2  2    -1.68642248E-02   # Re[N(2,2)]
   2  3    -7.06472551E-01   # Re[N(2,3)]
   2  4    -7.07474410E-01   # Re[N(2,4)]
   3  1     1.22258053E-01   # Re[N(3,1)]
   3  2     9.89502894E-01   # Re[N(3,2)]
   3  3     4.24411671E-02   # Re[N(3,3)]
   3  4    -6.43097115E-02   # Re[N(3,4)]
   4  1    -9.91384596E-01   # Re[N(4,1)]
   4  2     1.18075445E-01   # Re[N(4,2)]
   4  3     3.11546036E-02   # Re[N(4,3)]
   4  4    -4.73725947E-02   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -6.52600624E-02   # Re[U(1,1)]
   1  2     9.97868290E-01   # Re[U(1,2)]
   2  1     9.97868290E-01   # Re[U(2,1)]
   2  2     6.52600624E-02   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.88521167E-02   # Re[V(1,1)]
   1  2     9.95102135E-01   # Re[V(1,2)]
   2  1     9.95102135E-01   # Re[V(2,1)]
   2  2     9.88521167E-02   # Re[V(2,2)]
DECAY   2000011     1.53002210E-01   # ~e^-_R
#    BR                NDA      ID1      ID2
     6.78481910E-02    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
     2.91127992E-03    2     1000023        11   # BR(~e^-_R -> chi^0_2 e^-)
     1.63586739E-02    2     1000025        11   # BR(~e^-_R -> chi^0_3 e^-)
     9.12877540E-01    2     1000035        11   # BR(~e^-_R -> chi^0_4 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000011     4.67799697E-01   # ~e^-_L
#    BR                NDA      ID1      ID2
     2.47216310E-02    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     1.03306270E-03    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     3.13192607E-01    2     1000025        11   # BR(~e^-_L -> chi^0_3 e^-)
     4.74634639E-02    2     1000035        11   # BR(~e^-_L -> chi^0_4 e^-)
     6.84016161E-02    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
     5.45187619E-01    2    -1000037        12   # BR(~e^-_L -> chi^-_2 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.79630966E-01   # ~mu^-_R
#    BR                NDA      ID1      ID2
     7.14528951E-02    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     2.16506149E-03    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     6.88096876E-02    2     1000025        13   # BR(~mu^-_R -> chi^0_3 mu^-)
     7.30910570E-01    2     1000035        13   # BR(~mu^-_R -> chi^0_4 mu^-)
     2.96179033E-02    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
     9.70438824E-02    2    -1000037        14   # BR(~mu^-_R -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000013     4.43300484E-01   # ~mu^-_L
#    BR                NDA      ID1      ID2
     2.21576645E-02    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     2.81555124E-03    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     3.08264096E-01    2     1000025        13   # BR(~mu^-_L -> chi^0_3 mu^-)
     6.89910496E-02    2     1000035        13   # BR(~mu^-_L -> chi^0_4 mu^-)
     6.17770594E-02    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
     5.35994579E-01    2    -1000037        14   # BR(~mu^-_L -> chi^-_2 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     7.12678832E-01   # ~tau^-_1
#    BR                NDA      ID1      ID2
     2.48294278E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.21331407E-01    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     8.64254673E-02    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.13118052E-01    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     2.85135856E-01    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
     1.45694940E-01    2    -1000037        16   # BR(~tau^-_1 -> chi^-_2 nu_tau)
DECAY   2000015     5.11072174E-01   # ~tau^-_2
#    BR                NDA      ID1      ID2
     9.06419487E-02    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     2.24557400E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     1.71622755E-01    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     1.59556682E-01    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     5.66066891E-02    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
     2.97014481E-01    2    -1000037        16   # BR(~tau^-_2 -> chi^-_2 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     5.23413640E-01   # ~nu_e
#    BR                NDA      ID1      ID2
     8.36280848E-02    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     3.53994444E-03    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     2.05400758E-01    2     1000025        12   # BR(~nu_e -> chi^0_3 nu_e)
     9.61316274E-02    2     1000035        12   # BR(~nu_e -> chi^0_4 nu_e)
     1.39890201E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
     4.71409384E-01    2     1000037        11   # BR(~nu_e -> chi^+_2 e^-)
DECAY   1000014     5.24123576E-01   # ~nu_mu
#    BR                NDA      ID1      ID2
     8.35148098E-02    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     3.53514955E-03    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     2.05122563E-01    2     1000025        14   # BR(~nu_mu -> chi^0_3 nu_mu)
     9.60014279E-02    2     1000035        14   # BR(~nu_mu -> chi^0_4 nu_mu)
     1.41054639E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
     4.70771410E-01    2     1000037        13   # BR(~nu_mu -> chi^+_2 mu^-)
DECAY   1000016     7.24172434E-01   # ~nu_tau
#    BR                NDA      ID1      ID2
     6.04444818E-02    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     2.55859163E-03    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     1.48463599E-01    2     1000025        16   # BR(~nu_tau -> chi^0_3 nu_tau)
     6.94841346E-02    2     1000035        16   # BR(~nu_tau -> chi^0_4 nu_tau)
     3.78211822E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     3.40837371E-01    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000001     1.87038570E+02   # ~d_R
#    BR                NDA      ID1      ID2
     1.15006413E-04    2     1000025         1   # BR(~d_R -> chi^0_3 d)
     7.51288200E-03    2     1000035         1   # BR(~d_R -> chi^0_4 d)
     9.92348586E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     2.16177027E+02   # ~d_L
#    BR                NDA      ID1      ID2
     5.18626153E-04    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     4.30721045E-02    2     1000025         1   # BR(~d_L -> chi^0_3 d)
     4.30764056E-03    2     1000035         1   # BR(~d_L -> chi^0_4 d)
     5.40784694E-04    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     9.18009386E-02    2    -1000037         2   # BR(~d_L -> chi^-_2 u)
     8.59737763E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     1.87050006E+02   # ~s_R
#    BR                NDA      ID1      ID2
     1.31634303E-04    2     1000025         3   # BR(~s_R -> chi^0_3 s)
     7.51125530E-03    2     1000035         3   # BR(~s_R -> chi^0_4 s)
     9.92288392E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     2.16168453E+02   # ~s_L
#    BR                NDA      ID1      ID2
     5.18690726E-04    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     4.30594449E-02    2     1000025         3   # BR(~s_L -> chi^0_3 s)
     4.30883037E-03    2     1000035         3   # BR(~s_L -> chi^0_4 s)
     5.43306915E-04    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     9.17740196E-02    2    -1000037         4   # BR(~s_L -> chi^-_2 c)
     8.59771900E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     1.02646438E+02   # ~b_1
#    BR                NDA      ID1      ID2
     4.05077464E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     2.52765999E-03    2     1000023         5   # BR(~b_1 -> chi^0_2 b)
     8.79327301E-03    2     1000025         5   # BR(~b_1 -> chi^0_3 b)
     1.05298193E-03    2     1000035         5   # BR(~b_1 -> chi^0_4 b)
     1.57126369E-01    2    -1000024         6   # BR(~b_1 -> chi^-_1 t)
     1.88088687E-02    2    -1000037         6   # BR(~b_1 -> chi^-_2 t)
     6.05396617E-01    2     1000021         5   # BR(~b_1 -> ~g b)
     2.02243455E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     7.24774487E+01   # ~b_2
#    BR                NDA      ID1      ID2
     2.65297826E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     4.17136619E-03    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     2.81392202E-03    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     2.28619934E-03    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     5.26976476E-02    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     6.00408547E-03    2    -1000037         6   # BR(~b_2 -> chi^-_2 t)
     8.69070715E-01    2     1000021         5   # BR(~b_2 -> ~g b)
     6.03030861E-02    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   2000002     1.91310474E+02   # ~u_R
#    BR                NDA      ID1      ID2
     4.49535320E-04    2     1000025         2   # BR(~u_R -> chi^0_3 u)
     2.93758829E-02    2     1000035         2   # BR(~u_R -> chi^0_4 u)
     9.70082990E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     2.16162398E+02   # ~u_L
#    BR                NDA      ID1      ID2
     3.37174179E-04    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     4.73039469E-02    2     1000025         2   # BR(~u_L -> chi^0_3 u)
     2.26995731E-04    2     1000035         2   # BR(~u_L -> chi^0_4 u)
     1.24065024E-03    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     9.12704256E-02    2     1000037         1   # BR(~u_L -> chi^+_2 d)
     8.59606470E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     1.91312576E+02   # ~c_R
#    BR                NDA      ID1      ID2
     4.49638288E-04    2     1000025         4   # BR(~c_R -> chi^0_3 c)
     2.93754378E-02    2     1000035         4   # BR(~c_R -> chi^0_4 c)
     9.70073253E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     2.16163676E+02   # ~c_L
#    BR                NDA      ID1      ID2
     3.39091379E-04    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     4.73036193E-02    2     1000025         4   # BR(~c_L -> chi^0_3 c)
     2.27130255E-04    2     1000035         4   # BR(~c_L -> chi^0_4 c)
     1.24321284E-03    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     9.12698173E-02    2     1000037         3   # BR(~c_L -> chi^+_2 s)
     8.59600719E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     2.84268888E+01   # ~t_1
#    BR                NDA      ID1      ID2
     1.13514073E-01    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     1.32645552E-01    2     1000023         6   # BR(~t_1 -> chi^0_2 t)
     2.60764802E-01    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
     4.93068361E-01    2     1000021         6   # BR(~t_1 -> ~g t)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000006     1.09621632E+02   # ~t_2
#    BR                NDA      ID1      ID2
     9.05001550E-02    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     8.81451558E-02    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     1.18797470E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     4.38947340E-04    2     1000035         6   # BR(~t_2 -> chi^0_4 t)
     6.91773829E-03    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     2.31880800E-02    2     1000037         5   # BR(~t_2 -> chi^+_2 b)
     5.44884789E-01    2     1000021         6   # BR(~t_2 -> ~g t)
     1.23505205E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     1.10540182E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     1.60090710E-10   # chi^+_1
#    BR                NDA      ID1      ID2       ID3
     4.77924068E-01    3     1000022        -1         2   # BR(chi^+_1 -> chi^0_1 d_bar u)
     1.82863383E-01    3     1000022        -3         4   # BR(chi^+_1 -> chi^0_1 s_bar c)
     1.59315775E-01    3     1000022       -11        12   # BR(chi^+_1 -> chi^0_1 e^+ nu_e)
     1.58243659E-01    3     1000022       -13        14   # BR(chi^+_1 -> chi^0_1 mu^+ nu_mu)
     2.16531153E-02    3     1000022       -15        16   # BR(chi^+_1 -> chi^0_1 tau^+ nu_tau)
DECAY   1000037     1.46372489E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     2.44061656E-01    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     2.44458765E-01    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.51965109E-01    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     2.48384018E-01    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.53852088E-03    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     3.76832206E-03    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     3.70270180E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     3.92001985E-08   # chi^0_2
#    BR                NDA      ID1      ID2
     5.23019248E-03    2     1000022        22   # BR(chi^0_2 -> chi^0_1 photon)
#    BR                NDA      ID1      ID2       ID3
     1.08896701E-01    3     1000022         2        -2   # BR(chi^0_2 -> chi^0_1 u u_bar)
     9.28114115E-02    3     1000022         4        -4   # BR(chi^0_2 -> chi^0_1 c c_bar)
     1.41310948E-01    3     1000022         1        -1   # BR(chi^0_2 -> chi^0_1 d d_bar)
     1.41124630E-01    3     1000022         3        -3   # BR(chi^0_2 -> chi^0_1 s s_bar)
     3.24364424E-02    3     1000022        11       -11   # BR(chi^0_2 -> chi^0_1 e^- e^+)
     3.24098846E-02    3     1000022        13       -13   # BR(chi^0_2 -> chi^0_1 mu^- mu^+)
     2.38906107E-02    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
     1.94555889E-01    3     1000022        12       -12   # BR(chi^0_2 -> chi^0_1 nu_e nu_bar_e)
     4.33060084E-02    3     1000024         1        -2   # BR(chi^0_2 -> chi^+_1 d u_bar)
     4.33060084E-02    3    -1000024        -1         2   # BR(chi^0_2 -> chi^-_1 d_bar u)
     3.29947573E-02    3     1000024         3        -4   # BR(chi^0_2 -> chi^+_1 s c_bar)
     3.29947573E-02    3    -1000024        -3         4   # BR(chi^0_2 -> chi^-_1 s_bar c)
     1.44353856E-02    3     1000024        11       -12   # BR(chi^0_2 -> chi^+_1 e^- nu_bar_e)
     1.44353856E-02    3    -1000024       -11        12   # BR(chi^0_2 -> chi^-_1 e^+ nu_e)
     1.44073224E-02    3     1000024        13       -14   # BR(chi^0_2 -> chi^+_1 mu^- nu_bar_mu)
     1.44073224E-02    3    -1000024       -13        14   # BR(chi^0_2 -> chi^-_1 mu^+ nu_mu)
     8.52317150E-03    3     1000024        15       -16   # BR(chi^0_2 -> chi^+_1 tau^- nu_bar_tau)
     8.52317150E-03    3    -1000024       -15        16   # BR(chi^0_2 -> chi^-_1 tau^+ nu_tau)
DECAY   1000025     1.44634206E+01   # chi^0_3
#    BR                NDA      ID1      ID2
     2.82541272E-01    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     2.82541272E-01    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     8.85203477E-03    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     2.04067486E-01    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     2.02661462E-01    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     8.39855606E-03    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.03846931E-03    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     1.69636031E-03    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     4.02238501E-03    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     4.02238501E-03    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     5.12250607E+00   # chi^0_4
#    BR                NDA      ID1      ID2
     1.98542282E-02    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     1.98542282E-02    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     1.40558568E-01    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     1.40558568E-01    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     1.42436543E-02    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     3.23540721E-01    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     3.18551119E-01    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     1.34975562E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.85704581E-03    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     3.00550959E-03    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     2.21490395E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     2.21490395E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     7.30642877E-06   # ~g
#    BR                NDA      ID1      ID2
     5.61996269E-01    2     1000022        21   # BR(~g -> chi^0_1 g)
     4.17195947E-01    2     1000023        21   # BR(~g -> chi^0_2 g)
#    BR                NDA      ID1      ID2       ID3
     9.48787722E-03    3     1000022         5        -5   # BR(~g -> chi^0_1 b b_bar)
     7.48264265E-03    3     1000023         5        -5   # BR(~g -> chi^0_2 b b_bar)
     1.65996839E-03    3     1000024         5        -6   # BR(~g -> chi^+_1 b t_bar)
     1.65996839E-03    3    -1000024        -5         6   # BR(~g -> chi^-_1 b_bar t)
DECAY        25     3.62753106E-03   # h^0
#    BR                NDA      ID1      ID2
     2.75286211E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     7.76847422E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.53520471E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.62302142E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.47708917E-02    2           4        -4   # BR(h^0 -> c c_bar)
     8.74114526E-02    2          21        21   # BR(h^0 -> g g)
     3.14791322E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.05187528E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.05187528E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.05187528E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     3.68156349E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     3.68156349E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.05187528E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.05187528E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.05187528E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     3.68156349E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     3.68156349E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     8.32237674E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     8.32237674E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     7.13346578E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     4.75564385E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     3.32895070E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     3.32895070E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     3.32895070E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     3.32895070E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     3.32895070E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     4.14994920E+00   # H^0
#    BR                NDA      ID1      ID2
     3.54237137E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     1.00169867E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     2.87748343E-04    2           3        -3   # BR(H^0 -> s s_bar)
     6.53962843E-01    2           5        -5   # BR(H^0 -> b b_bar)
     2.31363707E-01    2           6        -6   # BR(H^0 -> t t_bar)
     4.59069641E-03    2          23        23   # BR(H^0 -> Z Z)
     8.72093194E-03    2          24       -24   # BR(H^0 -> W^+ W^-)
     4.59196225E-04    2          25        25   # BR(H^0 -> h^0 h^0)
DECAY        36     4.10102781E+00   # A^0
#    BR                NDA      ID1      ID2
     3.58492296E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.01373460E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     2.91204706E-04    2           3        -3   # BR(A^0 -> s s_bar)
     6.61817643E-01    2           5        -5   # BR(A^0 -> b b_bar)
     2.35901756E-01    2           6        -6   # BR(A^0 -> t t_bar)
     1.91932010E-04    2          25        23   # BR(A^0 -> h^0 Z)
DECAY        37     3.58995040E+00   # H^+
#    BR                NDA      ID1      ID2
     4.09905516E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.15911949E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.95068566E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.83168118E-01    2          -5         6   # BR(H^+ -> b_bar t)
     2.14289777E-04    2          25        24   # BR(H^+ -> h^0 W^+)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.01168417E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    9.99883158E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99883158E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    1.01168417E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.01168417E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    9.99883158E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999661E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    3.39339994E-07        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999661E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    3.39339994E-07        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15912185E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    6.47861746E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.24737980E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    3.39339994E-07        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999661E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    2.97025021E-04   # BR(b -> s gamma)
    2    1.59126786E-06   # BR(b -> s mu+ mu-)
    3    3.54200477E-05   # BR(b -> s nu nu)
    4    5.18794845E-11   # BR(Bd -> mu+ mu-)
    5    3.38905641E-09   # BR(Bs -> mu+ mu-)
    6    1.10461997E-04   # BR(B_u -> tau nu)
    7    9.98846771E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73607742E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.03065575E+01   # |Delta(M_Bs)| [ps^-1] 
   20    7.48501150E-15   # Delta(g-2)_electron/2
   21    4.32365366E-11   # Delta(g-2)_muon/2
   22   -5.19722823E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    5.44248137E-03   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87780140E-01   # C7
     0305 4422   00   1     3.14556249E-03   # C7
     0305 4322   00   1     5.71321692E-05   # C7'
     0305 6421   00   0    -9.48966443E-02   # C8
     0305 6421   00   1    -5.76150597E-03   # C8
     0305 6321   00   1    -1.27617018E-04   # C8'
 03051111 4133   00   0    -4.75945644E-01   # C9 e+e-
 03051111 4133   00   1     4.55989207E-04   # C9 e+e-
 03051111 4233   00   1    -1.51662488E-06   # C9' e+e-
 03051111 4137   00   0    -3.93297639E+00   # C10 e+e-
 03051111 4137   00   1     1.57779336E-03   # C10 e+e-
 03051111 4237   00   1     4.02225944E-05   # C10' e+e-
 03051313 4133   00   0    -4.75945644E-01   # C9 mu+mu-
 03051313 4133   00   1     4.55986052E-04   # C9 mu+mu-
 03051313 4233   00   1    -1.51663263E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.93297639E+00   # C10 mu+mu-
 03051313 4137   00   1     1.57779018E-03   # C10 mu+mu-
 03051313 4237   00   1     4.02226022E-05   # C10' mu+mu-
 03051212 4237   00   0     1.49046984E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -3.83179929E-04   # C11 nu_1 nu_1
 03051212 4137   00   1    -9.70591600E-06   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49046984E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -3.83179929E-04   # C11 nu_2 nu_2
 03051414 4137   00   1    -9.70591454E-06   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49046984E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -3.83179923E-04   # C11 nu_3 nu_3
 03051616 4137   00   1    -9.70550372E-06   # C11' nu_3 nu_3
