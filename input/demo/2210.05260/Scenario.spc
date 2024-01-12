# WARNING: 3-loop corrections to the Higgs mass require 2L corrections to yt and yb.  The latter will be enabled.
# WARNING: 3-loop corrections to the Higgs mass require 3L beta functions.  The latter will be enabled.
# SOFTSUSY4.1.9 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    4.1.9       # version number
Block MODSEL  # Select model
     1    0   # nonUniversal
Block SMINPUTS             # Standard Model inputs
     1    1.27908970e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.18400000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.18000000e+00   # mb(mb)
     6    1.73200000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    3.49000000e+01   # tanb, DRbar, Feynman gauge
Block EXTPAR               # non-universal SUSY breaking parameters
      0   -1.00000000e+00  # Set MX=MSUSY
      1    1.43790000e+03  # M_1(MX)
      2    2.73960000e+03  # M_2(MX)
      3    3.07950000e+03  # M_3(MX)
     11    2.96530000e+03  # At(MX)
     12    3.05070000e+03  # Ab(MX)
     13    2.88030000e+03  # Atau(MX)
     23   -1.88080000e+03  # mu(MX)
     26    3.74220000e+03  # mA(pole)
     31    4.03410000e+03  # meL(MX)
     32    4.03410000e+03  # mmuL(MX)
     33    3.62020000e+03  # mtauL(MX)
     34    4.07512000e+03  # meR(MX)
     35    4.07512000e+03  # mmuR(MX)
     36    2.60590000e+03  # mtauR(MX)
     41    1.77320000e+03  # mqL1(MX)
     42    1.77320000e+03  # mqL2(MX)
     43    2.17270000e+03  # mqL3(MX)
     44    1.81610000e+03  # muR(MX)
     45    1.81610000e+03  # mcR(MX)
     46    1.42430000e+03  # mtR(MX)
     47    1.92680000e+03  # mdR(MX)
     48    1.92680000e+03  # msR(MX)
     49    2.91300000e+03  # mbR(MX)
# SOFTSUSY-specific non SLHA information:
# mixing=0 Desired accuracy=1.00000000e-04 Achieved accuracy=2.90608940e-05
# 3-loop RGE corrections are on. 2-loop Yukawa/g3 thresholds are on
# 2-loop O(a_s^2) threshold corrections to mt are on
# 2-loop O(a_s^2) threshold corrections mb are on
# 2-loop O(a_s y^2) and O(y^4) threshold corrections to mb are off
# 2-loop O(y^4) threshold corrections to mtau are off
# 2-loop O(a_s^2) and O(a_s y^2) threshold corrections to a_s are off
# 2-loop SUSY QCD computation of squark/gluino pole masses are on.
# Expansion around gluino and squark pole masses.
# Matching scale=1.73200000e+02
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03680140e+01   # MW
        25     1.23982408e+02   # h0
        35     3.74294390e+03   # H0
        36     3.74220050e+03   # A0
        37     3.74348955e+03   # H+
   1000021     3.05926488e+03   # ~g
   1000022     1.43569019e+03   # ~neutralino(1)
   1000023     1.88440107e+03   # ~neutralino(2)
   1000024     1.88293570e+03   # ~chargino(1)
   1000025    -1.88556114e+03   # ~neutralino(3)
   1000035     2.73895917e+03   # ~neutralino(4)
   1000037     2.73912323e+03   # ~chargino(2)
   1000001     1.83556889e+03   # ~d_L
   1000002     1.83372496e+03   # ~u_L
   1000003     1.83556889e+03   # ~s_L
   1000004     1.83372496e+03   # ~c_L
   1000005     2.22697048e+03   # ~b_1
   1000006     1.44630670e+03   # ~t_1
   1000011     4.04148066e+03   # ~e_L
   1000012     4.04037896e+03   # ~nue_L
   1000013     4.04148066e+03   # ~mu_L
   1000014     4.04037896e+03   # ~numu_L
   1000015     2.61350236e+03   # ~stau_1
   1000016     3.62801218e+03   # ~nu_tau_L
   2000001     1.97937095e+03   # ~d_R
   2000002     1.86856742e+03   # ~u_R
   2000003     1.97937095e+03   # ~s_R
   2000004     1.86856742e+03   # ~c_R
   2000005     2.97219829e+03   # ~b_2
   2000006     2.24794943e+03   # ~t_2
   2000011     4.07688883e+03   # ~e_R
   2000013     4.07688883e+03   # ~mu_R
   2000015     3.62948876e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.90610179e-02       # alpha - evaluated at p^2=0
Block nmix                  # neutralino mixing matrix Q=1.75638820e+03
  1  1     9.97600513e-01   # N_{1,1}
  1  2    -2.37175003e-03   # N_{1,2}
  1  3    -5.53224485e-02   # N_{1,3}
  1  4    -4.15574146e-02   # N_{1,4}
  2  1     6.85269377e-02   # N_{2,1}
  2  2     6.14946415e-02   # N_{2,2}
  2  3     7.04338927e-01   # N_{2,3}
  2  4     7.03867278e-01   # N_{2,4}
  3  1    -9.70940058e-03   # N_{3,1}
  3  2     1.21856122e-02   # N_{3,2}
  3  3    -7.06843769e-01   # N_{3,3}
  3  4     7.07198080e-01   # N_{3,4}
  4  1    -1.73307923e-03   # N_{4,1}
  4  2     9.98030207e-01   # N_{4,2}
  4  3    -3.48997019e-02   # N_{4,3}
  4  4    -5.21029032e-02   # N_{4,4}
Block Umix                  # chargino U mixing matrix  at Q=1.75638820e+03
  1  1     4.93470863e-02   # U_{1,1}
  1  2     9.98781690e-01   # U_{1,2}
  2  1    -9.98781690e-01   # U_{2,1}
  2  2     4.93470863e-02   # U_{2,2}
Block Vmix                  # chargino V mixing matrix  at Q=1.75638820e+03
  1  1     7.36692048e-02   # V_{1,1}
  1  2    -9.97282732e-01   # V_{1,2}
  2  1    -9.97282732e-01   # V_{2,1}
  2  2    -7.36692048e-02   # V_{2,2}
Block stopmix               # stop mixing matrix at Q=1.75638820e+03
  1  1    -1.54811505e-01   # F_{11}
  1  2     9.87944026e-01   # F_{12}
  2  1     9.87944026e-01   # F_{21}
  2  2     1.54811505e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix at Q=1.75638820e+03
  1  1     9.98580557e-01   # F_{11}
  1  2    -5.32622872e-02   # F_{12}
  2  1     5.32622872e-02   # F_{21}
  2  2     9.98580557e-01   # F_{22}
Block staumix               # stau mixing matrix at Q=1.75638820e+03
  1  1    -1.80032554e-02   # F_{11}
  1  2     9.99837928e-01   # F_{12}
  2  1     9.99837928e-01   # F_{21}
  2  2     1.80032554e-02   # F_{22}
Block gauge Q= 1.75638820e+03  # SM gauge couplings
     1     3.62909347e-01   # g'(Q)MSSM DRbar
     2     6.34056327e-01   # g(Q)MSSM DRbar
     3     1.02265056e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.75638820e+03  
  3  3     8.29255643e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.75638820e+03  
  3  3     5.92038373e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.75638820e+03  
  3  3     3.35479589e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.75638820e+03 # Higgs mixing parameters
     1    -1.88080000e+03    # mu(Q)MSSM DRbar
     2     3.44254312e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44194264e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.60053487e+07    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.75638820e+03  # MSSM DRbar SUSY breaking parameters
     1     1.43790000e+03      # M_1(Q)
     2     2.73960000e+03      # M_2(Q)
     3     3.07950000e+03      # M_3(Q)
    21     1.05157429e+07      # mH1^2(Q)
    22    -3.46358297e+06      # mH2^2(Q)
    31     4.03410000e+03      # meL(Q)
    32     4.03410000e+03      # mmuL(Q)
    33     3.62020000e+03      # mtauL(Q)
    34     4.07512000e+03      # meR(Q)
    35     4.07512000e+03      # mmuR(Q)
    36     2.60590000e+03      # mtauR(Q)
    41     1.77319999e+03      # mqL1(Q)
    42     1.77319999e+03      # mqL2(Q)
    43     2.17269999e+03      # mqL3(Q)
    44     1.81609999e+03      # muR(Q)
    45     1.81609999e+03      # mcR(Q)
    46     1.42429999e+03      # mtR(Q)
    47     1.92679999e+03      # mdR(Q)
    48     1.92679999e+03      # msR(Q)
    49     2.91300000e+03      # mbR(Q)
Block au Q= 1.75638820e+03  
  1  1     4.08312986e-06      # Au(Q)MSSM DRbar
  2  2     4.08323051e-06      # Ac(Q)MSSM DRbar
  3  3     2.96530001e+03      # At(Q)MSSM DRbar
Block ad Q= 1.75638820e+03  
  1  1     2.47468652e-06      # Ad(Q)MSSM DRbar
  2  2     2.47498882e-06      # As(Q)MSSM DRbar
  3  3     3.05070000e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.75638820e+03  
  1  1     1.45038242e-06      # Ae(Q)MSSM DRbar
  2  2     1.45058397e-06      # Amu(Q)MSSM DRbar
  3  3     2.88030000e+03      # Atau(Q)MSSM DRbar
#     PDG         Width             
DECAY 1000021     2.63982196e+02   # Gluino decays
#     BR                NDA   PDG1         PDG2       Comments          
      4.93800765e-02    2      1         -1000001   # ~g -> d ~d_L*
      4.93800765e-02    2     -1          1000001   # ~g -> db ~d_L
      4.07489979e-02    2      1         -2000001   # ~g -> d ~d_R*
      4.07489979e-02    2     -1          2000001   # ~g -> db ~d_R
      4.94916971e-02    2      2         -1000002   # ~g -> u ~u_L*
      4.94916971e-02    2     -2          1000002   # ~g -> ub ~u_L
      4.73852838e-02    2      2         -2000002   # ~g -> u ~u_R*
      4.73852838e-02    2     -2          2000002   # ~g -> ub ~u_R
      4.93800765e-02    2      3         -1000003   # ~g -> s ~s_L*
      4.93800765e-02    2     -3          1000003   # ~g -> sb ~s_L
      4.07489979e-02    2      3         -2000003   # ~g -> s ~s_R*
      4.07489979e-02    2     -3          2000003   # ~g -> sb ~s_R
      4.94916938e-02    2      4         -1000004   # ~g -> c ~c_L*
      4.94916938e-02    2     -4          1000004   # ~g -> cb ~c_L
      4.73852804e-02    2      4         -2000004   # ~g -> c ~c_R*
      4.73852804e-02    2     -4          2000004   # ~g -> cb ~c_R
      2.66616347e-02    2      5         -1000005   # ~g -> b ~b_1*
      2.66616347e-02    2     -5          1000005   # ~g -> bb ~b_1
      3.76651948e-04    2      5         -2000005   # ~g -> b ~b_2*
      3.76651948e-04    2     -5          2000005   # ~g -> bb ~b_2
      7.57356462e-02    2      6         -1000006   # ~g -> t ~t_1*
      7.57356462e-02    2     -6          1000006   # ~g -> tb ~t_1
      2.32139634e-02    2      6         -2000006   # ~g -> t ~t_2*
      2.32139634e-02    2     -6          2000006   # ~g -> tb ~t_2
#
#     PDG         Width             
DECAY 1000001     4.10266291e-02   # SdownL decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      1          1000022   # ~d_L -> d ~chi_10
#
#     PDG         Width             
DECAY 2000001     2.57625001e-01   # SdownR decays
#     BR                NDA   PDG1         PDG2       Comments          
      9.99809868e-01    2      1          1000022   # ~d_R -> d ~chi_10
      1.86331619e-04    2      1          1000023   # ~d_R -> d ~chi_20
      3.80010656e-06    2      1          1000025   # ~d_R -> d ~chi_30
#
#     PDG         Width             
DECAY 1000002     3.88554843e-02   # SupL decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      2          1000022   # ~u_L -> u ~chi_10
#
#     PDG         Width             
DECAY 2000002     7.26796889e-01   # SupR decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      2          1000022   # ~u_R -> u ~chi_10
#
#     PDG         Width             
DECAY 1000003     4.10266286e-02   # SstrangeL decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      3          1000022   # ~s_L -> s ~chi_10
#
#     PDG         Width             
DECAY 2000003     2.57624999e-01   # SstrangeR decays
#     BR                NDA   PDG1         PDG2       Comments          
      9.99809868e-01    2      3          1000022   # ~s_R -> s ~chi_10
      1.86331591e-04    2      3          1000023   # ~s_R -> s ~chi_20
      3.80010596e-06    2      3          1000025   # ~s_R -> s ~chi_30
#
#     PDG         Width             
DECAY 1000004     3.88554297e-02   # ScharmL decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      4          1000022   # ~c_L -> c ~chi_10
#
#     PDG         Width             
DECAY 2000004     7.26796012e-01   # ScharmR decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.00000000e+00    2      4          1000022   # ~c_R -> c ~chi_10
#
#     PDG         Width             
DECAY 1000005     3.53044663e+01   # Sbottom1 decays
#     BR                NDA   PDG1         PDG2       Comments          
      7.36728777e-02    2      6         -1000024   # ~b_1 -> t ~chi_1-
      8.80013288e-01    2     -24         1000006   # ~b_1 -> W- ~t_1
      4.21702086e-03    2      5          1000022   # ~b_1 -> b ~chi_10
      2.11365688e-02    2      5          1000023   # ~b_1 -> b ~chi_20
      2.09602442e-02    2      5          1000025   # ~b_1 -> b ~chi_30
#
#     PDG         Width             
DECAY 2000005     2.65542485e+01   # Sbottom2 decays
#     BR                NDA   PDG1         PDG2       Comments          
      3.19926214e-01    2      6         -1000024   # ~b_2 -> t ~chi_1-
      1.49337135e-04    2      6         -1000037   # ~b_2 -> t ~chi_2-
      1.84178937e-02    2     -24         1000006   # ~b_2 -> W- ~t_1
      1.30380374e-01    2     -24         2000006   # ~b_2 -> W- ~t_2
      8.75459658e-02    2      25         1000005   # ~b_2 -> h ~b_1
      7.36207471e-02    2      1000005    23        # ~b_2 -> Z ~b_1
      3.84811253e-02    2      5          1000022   # ~b_2 -> b ~chi_10
      1.64841746e-01    2      5          1000023   # ~b_2 -> b ~chi_20
      1.66530341e-01    2      5          1000025   # ~b_2 -> b ~chi_30
      1.06256187e-04    2      5          1000035   # ~b_2 -> b ~chi_40
#
#     PDG         Width             
DECAY 1000006     0.00000000e+00   # Stop1 decays
#
#     PDG         Width             
DECAY 2000006     3.84041064e+01   # Stop2 decays
#     BR                NDA   PDG1         PDG2       Comments          
      4.54125938e-02    2      5          1000024   # ~t_2 -> b ~chi_1+
      4.42854305e-01    2      25         1000006   # ~t_2 -> h ~t_1
      4.35609093e-01    2      1000006    23        # ~t_2 -> Z ~t_1
      6.46716890e-03    2      6          1000022   # ~t_2 -> t ~chi_10
      2.91657421e-02    2      6          1000023   # ~t_2 -> t ~chi_20
      4.04910976e-02    2      6          1000025   # ~t_2 -> t ~chi_30
#
#     PDG         Width             
DECAY 1000011     1.82602703e+01   # SelectronL decays
#     BR                NDA   PDG1         PDG2       Comments          
      2.48842098e-03    2      12        -1000024   # ~e_L -> nu_e ~chi_1-
      5.16246540e-01    2      12        -1000037   # ~e_L -> nu_e ~chi_2-
      2.18593795e-01    2      11         1000022   # ~e_L -> e- ~chi_10
      5.39396662e-03    2      11         1000023   # ~e_L -> e- ~chi_20
      2.29145134e-05    2      11         1000025   # ~e_L -> e- ~chi_30
      2.57254362e-01    2      11         1000035   # ~e_L -> e- ~chi_40
#
#     PDG         Width             
DECAY 2000011     1.63784110e+01   # SelectronR decays
#     BR                NDA   PDG1         PDG2       Comments          
      9.96088752e-01    2      11         1000022   # ~e_R -> e- ~chi_10
      3.83001390e-03    2      11         1000023   # ~e_R -> e- ~chi_20
      7.99527451e-05    2      11         1000025   # ~e_R -> e- ~chi_30
      1.28114310e-06    2      11         1000035   # ~e_R -> e- ~chi_40
#
#     PDG         Width             
DECAY 1000013     1.82602702e+01   # SmuonL decays
#     BR                NDA   PDG1         PDG2       Comments          
      2.48842099e-03    2      14        -1000024   # ~mu_L -> nu_mu ~chi_1-
      5.16246541e-01    2      14        -1000037   # ~mu_L -> nu_mu ~chi_2-
      2.18593795e-01    2      13         1000022   # ~mu_L -> mu- ~chi_10
      5.39396662e-03    2      13         1000023   # ~mu_L -> mu- ~chi_20
      2.29145134e-05    2      13         1000025   # ~mu_L -> mu- ~chi_30
      2.57254362e-01    2      13         1000035   # ~mu_L -> mu- ~chi_40
#
#     PDG         Width             
DECAY 2000013     1.63784110e+01   # SmuonR decays
#     BR                NDA   PDG1         PDG2       Comments          
      9.96088752e-01    2      13         1000022   # ~mu_R -> mu- ~chi_10
      3.83001390e-03    2      13         1000023   # ~mu_R -> mu- ~chi_20
      7.99527451e-05    2      13         1000025   # ~mu_R -> mu- ~chi_30
      1.28114310e-06    2      13         1000035   # ~mu_R -> mu- ~chi_40
#
#     PDG         Width             
DECAY 1000012     1.82630298e+01   # Selectron sneutrino decays
#     BR                NDA   PDG1         PDG2       Comments          
      5.66364569e-03    2      11         1000024   # ~nu_eL -> e- ~chi_1+
      5.14035654e-01    2      11         1000037   # ~nu_eL -> e- ~chi_2+
      2.21917390e-01    2      12         1000022   # ~nu_eL -> nu_e ~chi_10
      2.35951216e-04    2      12         1000023   # ~nu_eL -> nu_e ~chi_20
      1.73534261e-04    2      12         1000025   # ~nu_eL -> nu_e ~chi_30
      2.57973825e-01    2      12         1000035   # ~nu_eL -> nu_e ~chi_40
#
#     PDG         Width             
DECAY 1000014     1.82630298e+01   # Smuon sneutrino decays
#     BR                NDA   PDG1         PDG2       Comments          
      5.66364569e-03    2      13         1000024   # ~nu_muL -> mu- ~chi_1+
      5.14035653e-01    2      13         1000037   # ~nu_muL -> mu- ~chi_2+
      2.21917391e-01    2      14         1000022   # ~nu_muL -> nu_mu ~chi_10
      2.35951216e-04    2      14         1000023   # ~nu_muL -> nu_mu ~chi_20
      1.73534261e-04    2      14         1000025   # ~nu_muL -> nu_mu ~chi_30
      2.57973825e-01    2      14         1000035   # ~nu_muL -> nu_mu ~chi_40
#
#     PDG         Width             
DECAY 1000015     4.55975967e+00   # Stau1 decays
#     BR                NDA   PDG1         PDG2       Comments          
      3.03496703e-01    2      16        -1000024   # ~tau_1- -> nu_tau ~chi_1-
      3.62995949e-01    2      15         1000022   # ~tau_1- -> tau- ~chi_10
      1.68704988e-01    2      15         1000023   # ~tau_1- -> tau- ~chi_20
      1.64802360e-01    2      15         1000025   # ~tau_1- -> tau- ~chi_30
#
#     PDG         Width             
DECAY 2000015     2.48777550e+01   # Stau2 decays
#     BR                NDA   PDG1         PDG2       Comments          
      9.12492512e-04    2      16        -1000024   # ~tau_2- -> nu_tau ~chi_1-
      2.15837847e-01    2      16        -1000037   # ~tau_2- -> nu_tau ~chi_2-
      2.19163537e-02    2      1000015    23        # ~tau_2- -> Z ~tau_1-
      2.33123707e-02    2      1000015    25        # ~tau_2- -> h ~tau_1-
      5.42419464e-01    2      15         1000022   # ~tau_2- -> tau- ~chi_10
      9.85386440e-02    2      15         1000023   # ~tau_2- -> tau- ~chi_20
      9.70538952e-02    2      15         1000025   # ~tau_2- -> tau- ~chi_30
      8.93291619e-06    2      15         1000035   # ~tau_2- -> tau- ~chi_40
#
#     PDG         Width             
DECAY 1000016     3.27523487e+03   # Stau sneutrino decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.03532645e-03    2      16         1000022   # ~nu_tauL -> nu_tau ~chi_10
      1.02887453e-06    2      16         1000023   # ~nu_tauL -> nu_tau ~chi_20
      8.17857827e-04    2      16         1000035   # ~nu_tauL -> nu_tau ~chi_40
      1.49940049e-03    2      15         1000024   # ~nu_tauL -> tau- ~chi_1+
      1.63078851e-03    2      15         1000037   # ~nu_tauL -> tau- ~chi_2+
      9.95014841e-01    2      24         1000015   # ~nu_tauL -> W+ ~tau_1-
#
#     PDG         Width             
DECAY 1000024     8.08732864e+00   # Chargino 1+ (lightest) decays
#     BR                NDA   PDG1         PDG2       Comments          
      3.88268638e-05    2     -1          1000002   # ~chi_1+ -> db ~u_L
      1.58089811e-05    2      2         -1000001   # ~chi_1+ -> u ~d_L*
      3.88268423e-05    2     -3          1000004   # ~chi_1+ -> sb ~s_L
      1.58077663e-05    2      4         -1000003   # ~chi_1+ -> c ~s_L*
      9.04616406e-01    2     -5          1000006   # ~chi_1+ -> bb ~t_1
      9.52743235e-02    2      24         1000022   # ~chi_1+ -> W+ ~chi_10
#
#     PDG         Width             
DECAY 1000037     6.26761268e+01   # Chargino 2+ (heaviest) decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.58824953e-01    2     -1          1000002   # ~chi_2+ -> db ~u_L
      1.58772549e-01    2      2         -1000001   # ~chi_2+ -> u ~d_L*
      1.58824953e-01    2     -3          1000004   # ~chi_2+ -> sb ~c_L
      1.58772527e-01    2      4         -1000003   # ~chi_2+ -> c ~s_L*
      1.77940440e-02    2     -5          1000006   # ~chi_2+ -> bb ~t_1
      5.27200106e-02    2     -5          2000006   # ~chi_2+ -> bb ~t_2
      5.30517677e-02    2      6         -1000005   # ~chi_2+ -> t ~b_1*
      1.40132646e-03    2      16        -1000015   # ~chi_2+ -> nu_tau ~tau_1+
      1.00482065e-05    2      24         1000022   # ~chi_2+ -> W+ ~chi_10
      5.59640114e-02    2      24         1000023   # ~chi_2+ -> W+ ~chi_20
      5.49146064e-02    2      24         1000025   # ~chi_2+ -> W+ ~chi_30
      5.83393892e-02    2      1000024    23        # ~chi_2+ -> Z ~chi_1+
      7.06098146e-02    2      1000024    25        # ~chi_2+ -> h ~chi_1+
#
#     PDG         Width             
DECAY 1000022     0.00000000e+00   # Neutralino1 decays
#
#     PDG         Width             
DECAY 1000023     6.62742821e+00   # Neutralino2 decays
#     BR                NDA   PDG1         PDG2       Comments          
      2.59129768e-05    2     -2          1000002   # ~chi_20 -> ub ~u_L
      1.32086818e-06    2     -2          2000002   # ~chi_20 -> ub ~u_R
      2.59129768e-05    2      2         -1000002   # ~chi_20 -> u ~u_L*
      1.32086818e-06    2      2         -2000002   # ~chi_20 -> u ~u_R*
      9.92892360e-06    2     -1          1000001   # ~chi_20 -> db ~d_L
      9.92892360e-06    2      1         -1000001   # ~chi_20 -> d ~d_L*
      2.59112404e-05    2     -4          1000004   # ~chi_20 -> cb ~c_L
      1.31994385e-06    2     -4          2000004   # ~chi_20 -> cb ~c_R
      2.59112404e-05    2      4         -1000004   # ~chi_20 -> c ~c_L*
      1.31994385e-06    2      4         -2000004   # ~chi_20 -> c ~c_R*
      9.92891801e-06    2     -3          1000003   # ~chi_20 -> sb ~s_L
      9.92891801e-06    2      3         -1000003   # ~chi_20 -> s ~s_L*
      4.43989675e-01    2      6         -1000006   # ~chi_20 -> t ~t_1*
      4.43989675e-01    2     -6          1000006   # ~chi_20 -> tb ~t_1
      2.49249967e-03    2      23         1000022   # ~chi_20 -> Z ~chi_10
      1.09379504e-01    2      25         1000022   # ~chi_20 -> h ~chi_10
#
#     PDG         Width             
DECAY 1000025     8.78711592e+00   # Neutralino3 decays
#     BR                NDA   PDG1         PDG2       Comments          
      4.58886165e-01    2      6         -1000006   # ~chi_30 -> t ~t_1*
      4.58886165e-01    2     -6          1000006   # ~chi_30 -> tb ~t_1
      8.06339318e-02    2      23         1000022   # ~chi_30 -> Z ~chi_10
      1.58923018e-03    2      25         1000022   # ~chi_30 -> h ~chi_10
#
#     PDG         Width             
DECAY 1000035     6.21211534e+01   # Neutralino4 decays
#     BR                NDA   PDG1         PDG2       Comments          
      8.01614483e-02    2     -2          1000002   # ~chi_40 -> ub ~u_L
      8.01614483e-02    2      2         -1000002   # ~chi_40 -> u ~u_L*
      8.00100346e-02    2     -1          1000001   # ~chi_40 -> db ~d_L
      8.00100346e-02    2      1         -1000001   # ~chi_40 -> d ~d_L*
      8.01614372e-02    2     -4          1000004   # ~chi_40 -> cb ~c_L
      8.01614372e-02    2      4         -1000004   # ~chi_40 -> c ~c_L*
      8.00100346e-02    2     -3          1000003   # ~chi_40 -> sb ~s_L
      8.00100346e-02    2      3         -1000003   # ~chi_40 -> s ~s_L*
      8.78394791e-03    2      6         -1000006   # ~chi_40 -> t ~t_1*
      2.33609470e-02    2      6         -2000006   # ~chi_40 -> t ~t_2*
      8.78394791e-03    2     -6          1000006   # ~chi_40 -> tb ~t_1
      2.33609470e-02    2     -6          2000006   # ~chi_40 -> tb ~t_2
      3.01437679e-02    2      5         -1000005   # ~chi_40 -> b ~b_1*
      3.01437679e-02    2     -5          1000005   # ~chi_40 -> bb ~b_1
      7.02573392e-04    2     -1000015    15        # ~chi_40 -> tau- ~tau_1+
      7.02573392e-04    2      1000015   -15        # ~chi_40 -> tau+ ~tau_1-
      5.64720263e-02    2      24        -1000024   # ~chi_40 -> W+ ~chi_1-
      5.64720263e-02    2     -24         1000024   # ~chi_40 -> W- ~chi_1+
      2.30164554e-06    2      23         1000022   # ~chi_40 -> Z ~chi_10
      2.30632631e-03    2      23         1000023   # ~chi_40 -> Z ~chi_20
      5.50187327e-02    2      23         1000025   # ~chi_40 -> Z ~chi_30
      1.84420561e-05    2      25         1000022   # ~chi_40 -> h ~chi_10
      6.07445027e-02    2      25         1000023   # ~chi_40 -> h ~chi_20
      2.29657112e-03    2      25         1000025   # ~chi_40 -> h ~chi_30
#
#     PDG         Width             
DECAY 25          5.24762751e-03   # light higgs decays
#     BR                NDA   PDG1         PDG2       Comments          
      4.21563738e-02    2      4         -4         # h -> c cb    
      2.93607304e-04    2      3         -3         # h -> s sb    
      7.26911710e-01    2      5         -5         # h -> b bb    
      1.61556539e-04    2      13        -13        # h -> mu- mu+ 
      4.48652583e-02    2      15        -15        # h-> tau- tau+
      1.94999807e-03    2      22         22        # h -> gamma gamma
      5.19102168e-02    2      21         21        # h -> gluon gluon
      1.19889610e-03    2      23         22        # h -> Z gamma 
      1.17680949e-01    2      24        -24        # h -> WW* -> W f f'bar
      1.28714306e-02    2      23         23        # h -> ZZ* -> Z f f'bar
#
#     PDG         Width             
DECAY 35          1.00643350e+02   # heavy higgs decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.85133943e-04    2      3         -3         # H -> s sb    
      8.30810308e-01    2      5         -5         # H -> b bb    
      1.40335857e-03    2      6         -6         # H -> t tb    
      3.00783565e-04    2      13        -13        # H -> mu- mu+ 
      8.36233718e-02    2      15        -15        # H -> tau- tau+
      8.36840788e-05    2      1000022    1000022   # H -> ~chi_10 ~chi_10
      2.52045649e-03    2      1000022    1000023   # H -> ~chi_10 ~chi_20
      1.03234107e-02    2      1000022    1000025   # H -> ~chi_10 ~chi_30
      1.25047849e-02    2      1000006   -1000006   # H -> ~t_1 ~t_1*
      2.91167754e-02    2      1000006   -2000006   # H -> ~t_1 ~t_2*
      2.91167754e-02    2      2000006   -1000006   # H -> ~t_2 ~t_1*
      8.75986604e-06    2      21         21        # H -> gluon gluon
#
#     PDG         Width             
DECAY 36          9.99423508e+01   # pseudoscalar higgs decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.86438166e-04    2      3         -3         # A -> s sb    
      8.36502844e-01    2      5         -5         # A -> b bb    
      1.41769454e-03    2      6         -6         # A -> t tb    
      3.02833622e-04    2      13        -13        # A -> mu- mu+ 
      8.41933948e-02    2      15        -15        # A -> tau- tau+
      1.87667113e-04    2      1000022    1000022   # A -> ~chi_10 ~chi_10
      1.04381073e-02    2      1000022    1000023   # A -> ~chi_10 ~chi_20
      2.51477422e-03    2      1000022    1000025   # A -> ~chi_10 ~chi_30
      3.21196729e-02    2      1000006    2000006   # A-> ~t_1 ~t_2*
      3.21196729e-02    2      2000006    1000006   # A-> ~t_2 ~t_1*
      1.56472610e-05    2      21         21        # A -> gluon gluon
#
#     PDG         Width             
DECAY 37          2.29478634e+02   # Charged higgs+ decays
#     BR                NDA   PDG1         PDG2       Comments          
      1.05138213e-04    2      4         -3         # H+ -> c sb   
      9.21082513e-01    2      6         -5         # H+ -> t bb   
      1.36405959e-04    2      14        -13        # H+ -> mu+ nu_mu
      4.00663044e-02    2      16        -15        # H+ -> tau+ nu_tau
      6.09901433e-03    2      1000022    1000024   # H+ -> ~chi_10 ~chi_1+
      3.25097333e-02    2      1000006   -1000005   # H+ -> ~t_1 b_1*
#
