# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 29.03.2023,  20:06
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
    3    2.77573705E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.30000000E+03  # M_1
    2    1.79840000E+03  # M_2
    3    1.74480000E+03  # M_3
   23    2.65080000E+03  # mu
   25    2.70000000E+01  # tan(beta)
   26    1.44150000E+03  # m_A, pole mass
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
   1    3.61469548E-01  # g'(Q)^DRbar
   2    6.33776151E-01  # g(Q)^DRbar
   3    1.04109822E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.68637401E-06   # Y_u(Q)^DRbar
  2  2     3.67723071E-03   # Y_c(Q)^DRbar
  3  3     8.45107042E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     3.11183195E-04   # Y_d(Q)^DRbar
  2  2     6.53484145E-03   # Y_s(Q)^DRbar
  3  3     3.55498111E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.37339154E-05   # Y_e(Q)^DRbar
  2  2     1.73133860E-02   # Y_mu(Q)^DRbar
  3  3     2.90847595E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.15122835E-05   # A_u(Q)^DRbar
  2  2     2.71943775E-08   # A_c(Q)^DRbar
  3  3    -4.84897155E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     3.21354115E-07   # A_d(Q)^DRbar
  2  2     1.53025901E-08   # A_s(Q)^DRbar
  3  3     2.81295447E-10   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.19425921E-06   # A_e(Q)^DRbar
  2  2     5.77587770E-09   # A_mu(Q)^DRbar
  3  3     3.43822682E-10   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    1.30000000E+03  # M_1
   2    1.79840000E+03  # M_2
   3    1.74480000E+03  # M_3
  21   -4.97997631E+06  # M^2_(H,d)
  22   -7.35664727E+06  # M^2_(H,u)
  31    1.81550000E+03  # M_(L,11)
  32    1.81550000E+03  # M_(L,22)
  33    1.81550000E+03  # M_(L,33)
  34    1.81550000E+03  # M_(E,11)
  35    1.81550000E+03  # M_(E,22)
  36    1.81550000E+03  # M_(E,33)
  41    2.18970000E+03  # M_(Q,11)
  42    2.18970000E+03  # M_(Q,22)
  43    2.09530000E+03  # M_(Q,33)
  44    2.18970000E+03  # M_(U,11)
  45    2.18970000E+03  # M_(U,22)
  46    1.38800000E+03  # M_(U,33)
  47    2.18970000E+03  # M_(D,11)
  48    2.18970000E+03  # M_(D,22)
  49    2.09530000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02873418E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.24605218E+02  # h0
        35     1.44063425E+03  # H0
        36     1.44150000E+03  # A0
        37     1.44358326E+03  # H+
   1000001     2.21674195E+03  # ~d_L
   2000001     2.21262154E+03  # ~d_R
   1000002     2.21560140E+03  # ~u_L
   2000002     2.21166413E+03  # ~u_R
   1000003     2.21684141E+03  # ~s_L
   2000003     2.21252016E+03  # ~s_R
   1000004     2.21560112E+03  # ~c_L
   2000004     2.21166409E+03  # ~c_R
   1000005     2.05521544E+03  # ~b_1
   2000005     2.13710556E+03  # ~b_2
   1000006     1.36296439E+03  # ~t_1
   2000006     2.12884979E+03  # ~t_2
   1000011     1.82083694E+03  # ~e_L-
   2000011     1.81923843E+03  # ~e_R-
   1000012     1.81885383E+03  # ~nu_eL
   1000013     1.81770625E+03  # ~mu_L-
   2000013     1.82234311E+03  # ~mu_R-
   1000014     1.81884626E+03  # ~nu_muL
   1000015     1.77981505E+03  # ~tau_1-
   2000015     1.85291440E+03  # ~tau_2-
   1000016     1.81670669E+03  # ~nu_tauL
   1000021     1.86959295E+03  # ~g
   1000022     1.30632766E+03  # ~chi_10
   1000023     1.82689224E+03  # ~chi_20
   1000025    -2.63613462E+03  # ~chi_30
   1000035     2.64031535E+03  # ~chi_40
   1000024     1.82708978E+03  # ~chi_1+
   1000037     2.64005203E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -3.81884321E-02   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    2.65080000E+03  # mu
   2    2.70000000E+01  # tan[beta](Q)
   3    2.43963565E+02  # v(Q)
   4    2.33413937E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     2.78115250E-01   # Re[R_st(1,1)]
   1  2     9.60547712E-01   # Re[R_st(1,2)]
   2  1    -9.60547712E-01   # Re[R_st(2,1)]
   2  2     2.78115250E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     8.58057766E-01   # Re[R_sb(1,1)]
   1  2     5.13553182E-01   # Re[R_sb(1,2)]
   2  1    -5.13553182E-01   # Re[R_sb(2,1)]
   2  2     8.58057766E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.90720771E-01   # Re[R_sta(1,1)]
   1  2     7.23121578E-01   # Re[R_sta(1,2)]
   2  1    -7.23121578E-01   # Re[R_sta(2,1)]
   2  2     6.90720771E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.99637429E-01   # Re[N(1,1)]
   1  2     2.24694986E-03   # Re[N(1,2)]
   1  3    -2.36390160E-02   # Re[N(1,3)]
   1  4     1.26948307E-02   # Re[N(1,4)]
   2  1    -4.12323191E-03   # Re[N(2,1)]
   2  2    -9.97492782E-01   # Re[N(2,2)]
   2  3     5.75146790E-02   # Re[N(2,3)]
   2  4    -4.10269562E-02   # Re[N(2,4)]
   3  1    -7.70673003E-03   # Re[N(3,1)]
   3  2     1.16975812E-02   # Re[N(3,2)]
   3  3     7.06824667E-01   # Re[N(3,3)]
   3  4     7.07250071E-01   # Re[N(3,4)]
   4  1    -2.54679379E-02   # Re[N(4,1)]
   4  2     6.97586430E-02   # Re[N(4,2)]
   4  3     7.04650373E-01   # Re[N(4,3)]
   4  4    -7.05657826E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.96556341E-01   # Re[U(1,1)]
   1  2     8.29183941E-02   # Re[U(1,2)]
   2  1     8.29183941E-02   # Re[U(2,1)]
   2  2     9.96556341E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.98219940E-01   # Re[V(1,1)]
   1  2     5.96401908E-02   # Re[V(1,2)]
   2  1     5.96401908E-02   # Re[V(2,1)]
   2  2     9.98219940E-01   # Re[V(2,2)]
DECAY   2000011     2.21741233E+00   # ~e^-_R
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000011     5.52625620E-01   # ~e^-_L
#    BR                NDA      ID1      ID2
     9.99999835E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000013     1.66282208E+00   # ~mu^-_L
#    BR                NDA      ID1      ID2
     9.99999990E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     1.10713819E+00   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99999925E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     1.27831374E+00   # ~tau^-_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
DECAY   2000015     1.46713869E+00   # ~tau^-_2
#    BR                NDA      ID1      ID2
     9.93958087E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     1.89604551E-03    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     3.72806315E-03    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     5.58128540E-01   # ~nu_e
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000014     5.58116316E-01   # ~nu_mu
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000016     5.54720060E-01   # ~nu_tau
#    BR                NDA      ID1      ID2
     9.99898453E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000001     1.09512227E+01   # ~d_R
#    BR                NDA      ID1      ID2
     4.94889249E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.50496750E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     1.34765042E+01   # ~d_L
#    BR                NDA      ID1      ID2
     1.03553342E-02    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     6.71891101E-02    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.34214091E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     7.88241464E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     1.10008347E+01   # ~s_R
#    BR                NDA      ID1      ID2
     4.84174431E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     1.92423016E-03    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     3.84415123E-03    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
     9.45814175E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     1.34267736E+01   # ~s_L
#    BR                NDA      ID1      ID2
     1.10886465E-02    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     6.58646517E-02    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.31566599E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     7.91480103E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     5.71364768E+01   # ~b_1
#    BR                NDA      ID1      ID2
     3.81913289E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     4.88708674E-03    2     1000023         5   # BR(~b_1 -> chi^0_2 b)
     5.46067407E-03    2    -1000024         6   # BR(~b_1 -> chi^-_1 t)
     6.26906180E-02    2     1000021         5   # BR(~b_1 -> ~g b)
     9.23142488E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     3.27157395E+01   # ~b_2
#    BR                NDA      ID1      ID2
     1.09781983E-02    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     4.23347347E-03    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     6.23023963E-03    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     2.04156129E-01    2     1000021         5   # BR(~b_2 -> ~g b)
     7.74401960E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   2000002     1.25246737E+01   # ~u_R
#    BR                NDA      ID1      ID2
     1.72858923E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     8.27140377E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     1.34063335E+01   # ~u_L
#    BR                NDA      ID1      ID2
     9.91107944E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     6.74275055E-02    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.34715833E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     7.87945582E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     1.25251219E+01   # ~c_R
#    BR                NDA      ID1      ID2
     1.72845592E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     8.27142608E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     1.34056820E+01   # ~c_L
#    BR                NDA      ID1      ID2
     9.91796713E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     6.74266342E-02    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.34715244E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     7.87940155E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     2.71207641E-05   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     8.27192005E+01   # ~t_2
#    BR                NDA      ID1      ID2
     2.28065624E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     4.60123719E-03    2     1000023         6   # BR(~t_2 -> chi^0_2 t)
     1.25315495E-02    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     3.66500519E-02    2     1000021         6   # BR(~t_2 -> ~g t)
     5.32499795E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     4.11436710E-01    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     5.62730666E-01   # chi^+_1
#    BR                NDA      ID1      ID2
     6.01482570E-04    2    -1000011        12   # BR(chi^+_1 -> ~e^+_L nu_e)
     4.47548445E-04    2    -1000013        14   # BR(chi^+_1 -> ~mu^+_L nu_mu)
     2.32105083E-04    2    -2000013        14   # BR(chi^+_1 -> ~mu^+_R nu_mu)
     1.73449751E-02    2    -1000015        16   # BR(chi^+_1 -> ~tau^+_1 nu_tau)
     1.04590332E-03    2     1000012       -11   # BR(chi^+_1 -> ~nu_e e^+)
     1.04780369E-03    2     1000014       -13   # BR(chi^+_1 -> ~nu_mu mu^+)
     1.65979616E-03    2     1000016       -15   # BR(chi^+_1 -> ~nu_tau tau^+)
     9.61406637E-01    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     1.58232646E-02    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     2.09930505E-04    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
     1.54720085E-04    3     1000022       -15        16   # BR(chi^+_1 -> chi^0_1 tau^+ nu_tau)
DECAY   1000037     5.42507776E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     3.67501494E-04    2    -1000011        12   # BR(chi^+_2 -> ~e^+_L nu_e)
     3.71942797E-04    2    -2000013        14   # BR(chi^+_2 -> ~mu^+_R nu_mu)
     4.32598698E-03    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     7.07167470E-03    2    -2000015        16   # BR(chi^+_2 -> ~tau^+_2 nu_tau)
     1.90880365E-04    2     1000012       -11   # BR(chi^+_2 -> ~nu_e e^+)
     2.30628318E-04    2     1000014       -13   # BR(chi^+_2 -> ~nu_mu mu^+)
     1.14559702E-02    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     1.81428996E-04    2     1000002        -1   # BR(chi^+_2 -> ~u_L d_bar)
     1.87018982E-04    2     1000004        -3   # BR(chi^+_2 -> ~c_L s_bar)
     5.07886191E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     3.33135809E-02    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     3.48965136E-04    2    -1000001         2   # BR(chi^+_2 -> ~d^*_L u)
     3.55327862E-04    2    -1000003         4   # BR(chi^+_2 -> ~s^*_L c)
     1.32234738E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     3.10404218E-02    2    -2000005         6   # BR(chi^+_2 -> ~b^*_2 t)
     3.14905427E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     7.96045702E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     7.78608530E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     7.76167418E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.32338136E-03    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     1.13261908E-03    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     1.30151987E-03    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     6.39313371E-01   # chi^0_2
#    BR                NDA      ID1      ID2
     2.49948815E-04    2     1000011       -11   # BR(chi^0_2 -> ~e^-_L e^+)
     2.49948815E-04    2    -1000011        11   # BR(chi^0_2 -> ~e^+_L e^-)
     1.89995594E-04    2     1000013       -13   # BR(chi^0_2 -> ~mu^-_L mu^+)
     1.89995594E-04    2    -1000013        13   # BR(chi^0_2 -> ~mu^+_L mu^-)
     7.62882777E-03    2     1000015       -15   # BR(chi^0_2 -> ~tau^-_1 tau^+)
     7.62882777E-03    2    -1000015        15   # BR(chi^0_2 -> ~tau^+_1 tau^-)
     4.35884661E-04    2     1000012       -12   # BR(chi^0_2 -> ~nu_e nu_bar_e)
     4.35884661E-04    2    -1000012        12   # BR(chi^0_2 -> ~nu^*_e nu_e)
     4.36704363E-04    2     1000014       -14   # BR(chi^0_2 -> ~nu_mu nu_bar_mu)
     4.36704363E-04    2    -1000014        14   # BR(chi^0_2 -> ~nu^*_mu nu_mu)
     6.99018640E-04    2     1000016       -16   # BR(chi^0_2 -> ~nu_tau nu_bar_tau)
     6.99018640E-04    2    -1000016        16   # BR(chi^0_2 -> ~nu^*_tau nu_tau)
     7.15607070E-02    2     1000006        -4   # BR(chi^0_2 -> ~t_1 c_bar)
     7.15607070E-02    2    -1000006         4   # BR(chi^0_2 -> ~t^*_1 c)
     4.12628187E-01    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     4.12628187E-01    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     2.60220065E-04    2     1000022        23   # BR(chi^0_2 -> chi^0_1 Z)
     1.16262931E-02    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.45057973E-04    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
DECAY   1000025     5.31685489E+01   # chi^0_3
#    BR                NDA      ID1      ID2
     5.94276052E-03    2     1000015       -15   # BR(chi^0_3 -> ~tau^-_1 tau^+)
     5.94276052E-03    2    -1000015        15   # BR(chi^0_3 -> ~tau^+_1 tau^-)
     5.53879637E-03    2     2000015       -15   # BR(chi^0_3 -> ~tau^-_2 tau^+)
     5.53879637E-03    2    -2000015        15   # BR(chi^0_3 -> ~tau^+_2 tau^-)
     2.56929632E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     2.56929632E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     7.22172827E-02    2     2000006        -6   # BR(chi^0_3 -> ~t_2 t_bar)
     7.22172827E-02    2    -2000006         6   # BR(chi^0_3 -> ~t^*_2 t)
     1.39107235E-02    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     1.39107235E-02    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     1.13476689E-02    2     2000005        -5   # BR(chi^0_3 -> ~b_2 b_bar)
     1.13476689E-02    2    -2000005         5   # BR(chi^0_3 -> ~b^*_2 b)
     7.71366157E-02    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     7.71366157E-02    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     3.08881117E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     7.33784106E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     2.84258454E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     2.14490054E-03    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     9.96398425E-04    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     7.94324728E-04    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     1.30510358E-03    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     1.30510358E-03    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     5.47157237E+01   # chi^0_4
#    BR                NDA      ID1      ID2
     1.38007088E-04    2     2000013       -13   # BR(chi^0_4 -> ~mu^-_R mu^+)
     1.38007088E-04    2    -2000013        13   # BR(chi^0_4 -> ~mu^+_R mu^-)
     4.96334689E-03    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     4.96334689E-03    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     6.18208160E-03    2     2000015       -15   # BR(chi^0_4 -> ~tau^-_2 tau^+)
     6.18208160E-03    2    -2000015        15   # BR(chi^0_4 -> ~tau^+_2 tau^-)
     1.89075813E-04    2     1000012       -12   # BR(chi^0_4 -> ~nu_e nu_bar_e)
     1.89075813E-04    2    -1000012        12   # BR(chi^0_4 -> ~nu^*_e nu_e)
     1.89078657E-04    2     1000014       -14   # BR(chi^0_4 -> ~nu_mu nu_bar_mu)
     1.89078657E-04    2    -1000014        14   # BR(chi^0_4 -> ~nu^*_mu nu_mu)
     1.89882527E-04    2     1000016       -16   # BR(chi^0_4 -> ~nu_tau nu_bar_tau)
     1.89882527E-04    2    -1000016        16   # BR(chi^0_4 -> ~nu^*_tau nu_tau)
     1.06667903E-04    2     1000002        -2   # BR(chi^0_4 -> ~u_L u_bar)
     1.06667903E-04    2    -1000002         2   # BR(chi^0_4 -> ~u^*_L u)
     1.07618210E-04    2     1000004        -4   # BR(chi^0_4 -> ~c_L c_bar)
     1.07618210E-04    2    -1000004         4   # BR(chi^0_4 -> ~c^*_L c)
     2.89590246E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     2.89590246E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     4.79640108E-02    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     4.79640108E-02    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     1.40188804E-04    2     1000001        -1   # BR(chi^0_4 -> ~d_L d_bar)
     1.40188804E-04    2    -1000001         1   # BR(chi^0_4 -> ~d^*_L d)
     1.46129813E-04    2     1000003        -3   # BR(chi^0_4 -> ~s_L s_bar)
     1.46129813E-04    2    -1000003         3   # BR(chi^0_4 -> ~s^*_L s)
     1.24316100E-02    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     1.24316100E-02    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     1.21288179E-02    2     2000005        -5   # BR(chi^0_4 -> ~b_2 b_bar)
     1.21288179E-02    2    -2000005         5   # BR(chi^0_4 -> ~b^*_2 b)
     6.94633567E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     6.94633567E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     2.75157190E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     2.06961043E-03    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     2.82369637E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     7.45965772E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     9.16033413E-04    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     5.19259972E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     1.32504428E-03    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     1.32504428E-03    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     6.89070354E+00   # ~g
#    BR                NDA      ID1      ID2
     1.36953078E-02    2     1000006        -4   # BR(~g -> ~t_1 c_bar)
     1.36953078E-02    2    -1000006         4   # BR(~g -> ~t^*_1 c)
     4.86230279E-01    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     4.86230279E-01    2    -1000006         6   # BR(~g -> ~t^*_1 t)
#    BR                NDA      ID1      ID2       ID3
DECAY        25     3.24593904E-03   # h^0
#    BR                NDA      ID1      ID2
     3.43162314E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     9.60741073E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     1.90532607E-04    2           3        -3   # BR(h^0 -> s s_bar)
     5.45574814E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.95235769E-02    2           4        -4   # BR(h^0 -> c c_bar)
     9.14957271E-02    2          21        21   # BR(h^0 -> g g)
     3.38343203E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.00749196E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.00749196E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.00749196E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     3.52622187E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     3.52622187E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.00749196E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.00749196E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.00749196E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     3.52622187E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     3.52622187E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     7.67045819E-04    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     7.67045819E-04    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     6.57467845E-04    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     4.38311897E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     3.06818328E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     3.06818328E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     3.06818328E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     3.06818328E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     3.06818328E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     1.49247086E+01   # H^0
#    BR                NDA      ID1      ID2
     5.74789093E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     1.62207712E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     2.80441175E-04    2           3        -3   # BR(H^0 -> s s_bar)
     8.29941754E-01    2           5        -5   # BR(H^0 -> b b_bar)
     6.51277388E-03    2           6        -6   # BR(H^0 -> t t_bar)
     1.46960831E-04    2          23        23   # BR(H^0 -> Z Z)
     2.79979379E-04    2          24       -24   # BR(H^0 -> W^+ W^-)
DECAY        36     1.49271842E+01   # A^0
#    BR                NDA      ID1      ID2
     5.75089695E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.62293528E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     2.80589672E-04    2           3        -3   # BR(A^0 -> s s_bar)
     8.30382253E-01    2           5        -5   # BR(A^0 -> b b_bar)
     6.36993563E-03    2           6        -6   # BR(A^0 -> t t_bar)
DECAY        37     1.30823907E+01   # H^+
#    BR                NDA      ID1      ID2
     6.57133418E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.85446727E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.80975736E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.13520293E-01    2          -5         6   # BR(H^+ -> b_bar t)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.06408271E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    7.28935917E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    7.29000000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99912095E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    1.45964707E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    1.37174211E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.06408271E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    7.28935917E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    7.29000000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99998635E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    1.36496217E-06        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99998635E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    1.36496217E-06        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.12656860E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    8.84303028E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.17061258E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    1.36496217E-06        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99998635E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    2.93230971E-04   # BR(b -> s gamma)
    2    1.59171127E-06   # BR(b -> s mu+ mu-)
    3    3.54345450E-05   # BR(b -> s nu nu)
    4    5.43068282E-11   # BR(Bd -> mu+ mu-)
    5    3.55394540E-09   # BR(Bs -> mu+ mu-)
    6    1.09053950E-04   # BR(B_u -> tau nu)
    7    9.86114582E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.74169821E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.02833586E+01   # |Delta(M_Bs)| [ps^-1] 
   20    2.94584632E-14   # Delta(g-2)_electron/2
   21    8.68297435E-11   # Delta(g-2)_muon/2
   22   -8.96318982E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    1.33605842E-02   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87796016E-01   # C7
     0305 4422   00   1     7.00968114E-03   # C7
     0305 4322   00   1     1.27223987E-04   # C7'
     0305 6421   00   0    -9.49017158E-02   # C8
     0305 6421   00   1    -1.10211814E-03   # C8
     0305 6321   00   1    -2.27304277E-05   # C8'
 03051111 4133   00   0    -4.76639521E-01   # C9 e+e-
 03051111 4133   00   1     2.02900758E-04   # C9 e+e-
 03051111 4233   00   1    -1.47900739E-05   # C9' e+e-
 03051111 4137   00   0    -3.93307115E+00   # C10 e+e-
 03051111 4137   00   1     1.26259133E-03   # C10 e+e-
 03051111 4237   00   1     4.33141751E-04   # C10' e+e-
 03051313 4133   00   0    -4.76639521E-01   # C9 mu+mu-
 03051313 4133   00   1     2.02845706E-04   # C9 mu+mu-
 03051313 4233   00   1    -1.47909983E-05   # C9' mu+mu-
 03051313 4137   00   0    -3.93307115E+00   # C10 mu+mu-
 03051313 4137   00   1     1.26253623E-03   # C10 mu+mu-
 03051313 4237   00   1     4.33142675E-04   # C10' mu+mu-
 03051212 4237   00   0     1.49069574E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -3.07964002E-04   # C11 nu_1 nu_1
 03051212 4137   00   1    -1.04538573E-04   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49069574E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -3.07964002E-04   # C11 nu_2 nu_2
 03051414 4137   00   1    -1.04538357E-04   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49069574E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -3.07963944E-04   # C11 nu_3 nu_3
 03051616 4137   00   1    -1.04477392E-04   # C11' nu_3 nu_3
