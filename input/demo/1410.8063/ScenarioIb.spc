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
    3    1.68333882E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    9.89000000E+02  # M_1
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
   1    3.61007154E-01  # g'(Q)^DRbar
   2    6.32141158E-01  # g(Q)^DRbar
   3    1.01702769E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.56907445E-06   # Y_u(Q)^DRbar
  2  2     3.62757355E-03   # Y_c(Q)^DRbar
  3  3     8.40794280E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.11152914E-04   # Y_d(Q)^DRbar
  2  2     4.43420768E-03   # Y_s(Q)^DRbar
  3  3     2.13887764E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     4.92468352E-05   # Y_e(Q)^DRbar
  2  2     1.01826293E-02   # Y_mu(Q)^DRbar
  3  3     1.71167070E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.16698718E-05   # A_u(Q)^DRbar
  2  2     2.75666361E-08   # A_c(Q)^DRbar
  3  3    -4.38169013E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     4.73590433E-07   # A_d(Q)^DRbar
  2  2     2.25519433E-08   # A_s(Q)^DRbar
  3  3     4.67534927E-10   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.03058734E-06   # A_e(Q)^DRbar
  2  2     9.82064619E-09   # A_mu(Q)^DRbar
  3  3     5.84224525E-10   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    9.89000000E+02  # M_1
   2    1.94410000E+03  # M_2
   3    5.83240000E+03  # M_3
  21   -3.38524631E+06  # M^2_(H,d)
  22   -7.38391136E+06  # M^2_(H,u)
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
  46    2.15900000E+03  # M_(U,33)
  47    3.05430000E+03  # M_(D,11)
  48    3.05430000E+03  # M_(D,22)
  49    2.14370000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02835287E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.29381007E+02  # h0
        35     1.91748035E+03  # H0
        36     1.91790000E+03  # A0
        37     1.91954610E+03  # H+
   1000001     2.57932711E+03  # ~d_L
   2000001     2.58185898E+03  # ~d_R
   1000002     2.57817085E+03  # ~u_L
   2000002     2.58153028E+03  # ~u_R
   1000003     2.57926897E+03  # ~s_L
   2000003     2.58191819E+03  # ~s_R
   1000004     2.57817097E+03  # ~c_L
   2000004     2.58153143E+03  # ~c_R
   1000005     1.26575692E+03  # ~b_1
   2000005     1.34012671E+03  # ~b_2
   1000006     1.07964295E+03  # ~t_1
   2000006     1.54344338E+03  # ~t_2
   1000011     2.25248403E+03  # ~e_L-
   2000011     2.24991924E+03  # ~e_R-
   1000012     2.25077935E+03  # ~nu_eL
   1000013     2.25285221E+03  # ~mu_L-
   2000013     2.24954749E+03  # ~mu_R-
   1000014     2.25077854E+03  # ~nu_muL
   1000015     2.23323428E+03  # ~tau_1-
   2000015     2.26832341E+03  # ~tau_2-
   1000016     2.25054891E+03  # ~nu_tauL
   1000021     5.28840804E+03  # ~g
   1000022     9.99567509E+02  # ~chi_10
   1000023     1.98682558E+03  # ~chi_20
   1000025    -2.64636466E+03  # ~chi_30
   1000035     2.65141734E+03  # ~chi_40
   1000024     1.98690429E+03  # ~chi_1+
   1000037     2.65142951E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -6.17709929E-02   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    2.65310000E+03  # mu
   2    1.63000000E+01  # tan[beta](Q)
   3    2.45311722E+02  # v(Q)
   4    2.21598273E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     7.37079483E-01   # Re[R_st(1,1)]
   1  2     6.75806064E-01   # Re[R_st(1,2)]
   2  1    -6.75806064E-01   # Re[R_st(2,1)]
   2  2     7.37079483E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     6.55472947E-01   # Re[R_sb(1,1)]
   1  2     7.55218654E-01   # Re[R_sb(1,2)]
   2  1    -7.55218654E-01   # Re[R_sb(2,1)]
   2  2     6.55472947E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.80057221E-01   # Re[R_sta(1,1)]
   1  2     7.33159039E-01   # Re[R_sta(1,2)]
   2  1    -7.33159039E-01   # Re[R_sta(2,1)]
   2  2     6.80057221E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99757138E-01   # Re[N(1,1)]
   1  2    -8.38104813E-04   # Re[N(1,2)]
   1  3     2.01806272E-02   # Re[N(1,3)]
   1  4    -8.81508169E-03   # Re[N(1,4)]
   2  1    -2.67183101E-03   # Re[N(2,1)]
   2  2    -9.96297401E-01   # Re[N(2,2)]
   2  3     6.80870549E-02   # Re[N(2,3)]
   2  4    -5.24261739E-02   # Re[N(2,4)]
   3  1     8.02392943E-03   # Re[N(3,1)]
   3  2    -1.11163474E-02   # Re[N(3,2)]
   3  3    -7.06881614E-01   # Re[N(3,3)]
   3  4    -7.07199001E-01   # Re[N(3,4)]
   4  1    -2.03505203E-02   # Re[N(4,1)]
   4  2     8.52479550E-02   # Re[N(4,2)]
   4  3     7.03757969E-01   # Re[N(4,3)]
   4  4    -7.05013024E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.95113555E-01   # Re[U(1,1)]
   1  2     9.87370930E-02   # Re[U(1,2)]
   2  1     9.87370930E-02   # Re[U(2,1)]
   2  2     9.95113555E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.97063633E-01   # Re[V(1,1)]
   1  2     7.65774905E-02   # Re[V(1,2)]
   2  1     7.65774905E-02   # Re[V(2,1)]
   2  2     9.97063633E-01   # Re[V(2,2)]
DECAY   2000011     7.51230478E+00   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99998791E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000011     3.18924796E+00   # ~e^-_L
#    BR                NDA      ID1      ID2
     5.88485228E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     1.37720287E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     2.73794485E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000013     7.02428230E+00   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.79272026E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
     6.93611046E-03    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
     1.37918633E-02    2    -1000024        14   # BR(~mu^-_R -> chi^-_1 nu_mu)
DECAY   1000013     3.67722355E+00   # ~mu^-_L
#    BR                NDA      ID1      ID2
     6.82695836E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     1.06193629E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     2.11110536E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000015     5.42233173E+00   # ~tau^-_1
#    BR                NDA      ID1      ID2
     8.96375871E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     3.46487099E-02    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
     6.89754189E-02    2    -1000024        16   # BR(~tau^-_1 -> chi^-_1 nu_tau)
DECAY   2000015     5.26819015E+00   # ~tau^-_2
#    BR                NDA      ID1      ID2
     8.58877507E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     4.72873163E-02    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     9.38331308E-02    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     3.18329744E+00   # ~nu_e
#    BR                NDA      ID1      ID2
     5.92161270E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     1.35580016E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     2.72258714E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000014     3.18329008E+00   # ~nu_mu
#    BR                NDA      ID1      ID2
     5.92162214E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     1.35579589E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     2.72258197E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000016     3.18121338E+00   # ~nu_tau
#    BR                NDA      ID1      ID2
     5.92428916E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     1.35458925E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     2.72111807E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000001     5.30864065E+00   # ~d_L
#    BR                NDA      ID1      ID2
     5.09732530E-02    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     3.16677209E-01    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     6.32349538E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
DECAY   2000001     1.07477622E+00   # ~d_R
#    BR                NDA      ID1      ID2
     9.99746811E-01    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     1.67595443E-04    2    -1000024         2   # BR(~d_R -> chi^-_1 u)
DECAY   1000003     5.21549652E+00   # ~s_L
#    BR                NDA      ID1      ID2
     5.53275242E-02    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     3.15224077E-01    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     6.29448399E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
DECAY   2000003     1.16790817E+00   # ~s_R
#    BR                NDA      ID1      ID2
     9.04645518E-01    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     3.18247072E-02    2     1000023         3   # BR(~s_R -> chi^0_2 s)
     6.35297753E-02    2    -1000024         4   # BR(~s_R -> chi^-_1 c)
DECAY   1000005     4.40099131E+00   # ~b_1
#    BR                NDA      ID1      ID2
     1.66453408E-02    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     9.83354659E-01    2     1000006       -24   # BR(~b_1 -> ~t_1 W^-)
DECAY   2000005     1.71719983E+01   # ~b_2
#    BR                NDA      ID1      ID2
     4.78827816E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     9.95211722E-01    2     1000006       -24   # BR(~b_2 -> ~t_1 W^-)
DECAY   1000002     5.30505770E+00   # ~u_L
#    BR                NDA      ID1      ID2
     5.00705987E-02    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     3.16582960E-01    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     6.33346442E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
DECAY   2000002     4.29726550E+00   # ~u_R
#    BR                NDA      ID1      ID2
     9.99998356E-01    2     1000022         2   # BR(~u_R -> chi^0_1 u)
DECAY   1000004     5.30482149E+00   # ~c_L
#    BR                NDA      ID1      ID2
     5.01862670E-02    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     3.16543497E-01    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     6.33270236E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
DECAY   2000004     4.29750395E+00   # ~c_R
#    BR                NDA      ID1      ID2
     9.99802534E-01    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     1.30635713E-04    2     1000024         3   # BR(~c_R -> chi^+_1 s)
DECAY   1000006     3.67823679E-09   # ~t_1
#    BR                NDA      ID1      ID2
     1.00000000E+00    2     1000022         4   # BR(~t_1 -> chi^0_1 c)
DECAY   2000006     5.76158208E+01   # ~t_2
#    BR                NDA      ID1      ID2
     9.38306169E-03    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     2.37591769E-01    2     1000005        24   # BR(~t_2 -> ~b_1 W^+)
     1.18366934E-01    2     2000005        24   # BR(~t_2 -> ~b_2 W^+)
     6.28968570E-01    2     1000006        23   # BR(~t_2 -> ~t_1 Z)
     5.68966565E-03    2     1000006        25   # BR(~t_2 -> ~t_1 h^0)
DECAY   1000024     1.69280076E+01   # chi^+_1
#    BR                NDA      ID1      ID2
     4.49658733E-01    2     1000006        -5   # BR(chi^+_1 -> ~t_1 b_bar)
     7.90858392E-02    2     2000006        -5   # BR(chi^+_1 -> ~t_2 b_bar)
     2.37546523E-01    2    -1000005         6   # BR(chi^+_1 -> ~b^*_1 t)
     2.32743492E-01    2    -2000005         6   # BR(chi^+_1 -> ~b^*_2 t)
     7.98115072E-04    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
     1.13121434E-04    3     1000022        -5         6   # BR(chi^+_1 -> chi^0_1 b_bar t)
DECAY   1000037     7.88695307E+01   # chi^+_2
#    BR                NDA      ID1      ID2
     1.00890854E-04    2    -1000011        12   # BR(chi^+_2 -> ~e^+_L nu_e)
     1.00041653E-04    2    -1000013        14   # BR(chi^+_2 -> ~mu^+_L nu_mu)
     1.91858974E-04    2    -1000015        16   # BR(chi^+_2 -> ~tau^+_1 nu_tau)
     6.27672036E-04    2    -2000015        16   # BR(chi^+_2 -> ~tau^+_2 nu_tau)
     8.17225404E-04    2     1000016       -15   # BR(chi^+_2 -> ~nu_tau tau^+)
     2.14367322E-01    2     1000006        -5   # BR(chi^+_2 -> ~t_1 b_bar)
     1.94303225E-01    2     2000006        -5   # BR(chi^+_2 -> ~t_2 b_bar)
     2.01755957E-01    2    -1000005         6   # BR(chi^+_2 -> ~b^*_1 t)
     2.18458348E-01    2    -2000005         6   # BR(chi^+_2 -> ~b^*_2 t)
     2.17757570E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     4.92570055E-02    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     4.75741002E-02    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     4.85830700E-02    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     4.99135407E-04    3     1000024         6        -6   # BR(chi^+_2 -> chi^+_1 t t_bar)
     7.79955566E-04    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     6.41506443E-04    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     1.68928404E+01   # chi^0_2
#    BR                NDA      ID1      ID2
     2.32957701E-01    2     1000006        -6   # BR(chi^0_2 -> ~t_1 t_bar)
     2.32957701E-01    2    -1000006         6   # BR(chi^0_2 -> ~t^*_1 t)
     4.14428826E-02    2     2000006        -6   # BR(chi^0_2 -> ~t_2 t_bar)
     4.14428826E-02    2    -2000006         6   # BR(chi^0_2 -> ~t^*_2 t)
     1.13715016E-01    2     1000005        -5   # BR(chi^0_2 -> ~b_1 b_bar)
     1.13715016E-01    2    -1000005         5   # BR(chi^0_2 -> ~b^*_1 b)
     1.11414174E-01    2     2000005        -5   # BR(chi^0_2 -> ~b_2 b_bar)
     1.11414174E-01    2    -2000005         5   # BR(chi^0_2 -> ~b^*_2 b)
     6.98191470E-04    2     1000022        25   # BR(chi^0_2 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.35293213E-04    3     1000022         6        -6   # BR(chi^0_2 -> chi^0_1 t t_bar)
DECAY   1000025     7.86223152E+01   # chi^0_3
#    BR                NDA      ID1      ID2
     3.81234125E-04    2     1000015       -15   # BR(chi^0_3 -> ~tau^-_1 tau^+)
     3.81234125E-04    2    -1000015        15   # BR(chi^0_3 -> ~tau^+_1 tau^-)
     3.67074320E-04    2     2000015       -15   # BR(chi^0_3 -> ~tau^-_2 tau^+)
     3.67074320E-04    2    -2000015        15   # BR(chi^0_3 -> ~tau^+_2 tau^-)
     2.10326149E-01    2     1000006        -6   # BR(chi^0_3 -> ~t_1 t_bar)
     2.10326149E-01    2    -1000006         6   # BR(chi^0_3 -> ~t^*_1 t)
     1.81580801E-01    2     2000006        -6   # BR(chi^0_3 -> ~t_2 t_bar)
     1.81580801E-01    2    -2000006         6   # BR(chi^0_3 -> ~t^*_2 t)
     1.29846548E-02    2     1000005        -5   # BR(chi^0_3 -> ~b_1 b_bar)
     1.29846548E-02    2    -1000005         5   # BR(chi^0_3 -> ~b^*_1 b)
     1.33452240E-02    2     2000005        -5   # BR(chi^0_3 -> ~b_2 b_bar)
     1.33452240E-02    2    -2000005         5   # BR(chi^0_3 -> ~b^*_2 b)
     4.61449822E-02    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     4.61449822E-02    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     1.97126730E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     4.38012538E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     3.15536124E-03    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     7.43504232E-04    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     6.92952696E-04    3     1000022         6        -6   # BR(chi^0_3 -> chi^0_1 t t_bar)
     3.36185668E-04    3     1000023         6        -6   # BR(chi^0_3 -> chi^0_2 t t_bar)
     6.36426248E-04    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     6.36426248E-04    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
DECAY   1000035     7.90637091E+01   # chi^0_4
#    BR                NDA      ID1      ID2
     2.79200007E-04    2     1000015       -15   # BR(chi^0_4 -> ~tau^-_1 tau^+)
     2.79200007E-04    2    -1000015        15   # BR(chi^0_4 -> ~tau^+_1 tau^-)
     4.85345445E-04    2     2000015       -15   # BR(chi^0_4 -> ~tau^-_2 tau^+)
     4.85345445E-04    2    -2000015        15   # BR(chi^0_4 -> ~tau^+_2 tau^-)
     2.61941497E-01    2     1000006        -6   # BR(chi^0_4 -> ~t_1 t_bar)
     2.61941497E-01    2    -1000006         6   # BR(chi^0_4 -> ~t^*_1 t)
     1.32082117E-01    2     2000006        -6   # BR(chi^0_4 -> ~t_2 t_bar)
     1.32082117E-01    2    -2000006         6   # BR(chi^0_4 -> ~t^*_2 t)
     1.01198636E-02    2     1000005        -5   # BR(chi^0_4 -> ~b_1 b_bar)
     1.01198636E-02    2    -1000005         5   # BR(chi^0_4 -> ~b^*_1 b)
     1.65828496E-02    2     2000005        -5   # BR(chi^0_4 -> ~b_2 b_bar)
     1.65828496E-02    2    -2000005         5   # BR(chi^0_4 -> ~b^*_2 b)
     4.16017075E-02    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     4.16017075E-02    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     3.02087895E-03    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     7.76736007E-04    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     1.95970536E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     4.76870835E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     7.85169784E-04    3     1000022         6        -6   # BR(chi^0_4 -> chi^0_1 t t_bar)
     1.89771528E-04    3     1000023         6        -6   # BR(chi^0_4 -> chi^0_2 t t_bar)
     6.46608198E-04    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     6.46608198E-04    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
DECAY   1000021     8.88649899E+02   # ~g
#    BR                NDA      ID1      ID2
     3.55833261E-02    2     1000002        -2   # BR(~g -> ~u_L u_bar)
     3.55833261E-02    2    -1000002         2   # BR(~g -> ~u^*_L u)
     3.55254904E-02    2     2000002        -2   # BR(~g -> ~u_R u_bar)
     3.55254904E-02    2    -2000002         2   # BR(~g -> ~u^*_R u)
     3.55838796E-02    2     1000004        -4   # BR(~g -> ~c_L c_bar)
     3.55838796E-02    2    -1000004         4   # BR(~g -> ~c^*_L c)
     3.55249121E-02    2     2000004        -4   # BR(~g -> ~c_R c_bar)
     3.55249121E-02    2    -2000004         4   # BR(~g -> ~c^*_R c)
     5.23997223E-02    2     1000006        -6   # BR(~g -> ~t_1 t_bar)
     5.23997223E-02    2    -1000006         6   # BR(~g -> ~t^*_1 t)
     5.48837018E-02    2     2000006        -6   # BR(~g -> ~t_2 t_bar)
     5.48837018E-02    2    -2000006         6   # BR(~g -> ~t^*_2 t)
     3.55634220E-02    2     1000001        -1   # BR(~g -> ~d_L d_bar)
     3.55634220E-02    2    -1000001         1   # BR(~g -> ~d^*_L d)
     3.55198312E-02    2     2000001        -1   # BR(~g -> ~d_R d_bar)
     3.55198312E-02    2    -2000001         1   # BR(~g -> ~d^*_R d)
     3.55638785E-02    2     1000003        -3   # BR(~g -> ~s_L s_bar)
     3.55638785E-02    2    -1000003         3   # BR(~g -> ~s^*_L s)
     3.55193553E-02    2     2000003        -3   # BR(~g -> ~s_R s_bar)
     3.55193553E-02    2    -2000003         3   # BR(~g -> ~s^*_R s)
     5.43236975E-02    2     1000005        -5   # BR(~g -> ~b_1 b_bar)
     5.43236975E-02    2    -1000005         5   # BR(~g -> ~b^*_1 b)
     5.37088396E-02    2     2000005        -5   # BR(~g -> ~b_2 b_bar)
     5.37088396E-02    2    -2000005         5   # BR(~g -> ~b^*_2 b)
#    BR                NDA      ID1      ID2       ID3
     1.00008674E-04    3     1000024         5        -6   # BR(~g -> chi^+_1 b t_bar)
     1.00008674E-04    3    -1000024        -5         6   # BR(~g -> chi^-_1 b_bar t)
DECAY        25     3.75036220E-03   # h^0
#    BR                NDA      ID1      ID2
     2.81048692E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     7.91419135E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     2.04982364E-04    2           3        -3   # BR(h^0 -> s s_bar)
     4.65011907E-01    2           5        -5   # BR(h^0 -> b b_bar)
     3.40643592E-02    2           4        -4   # BR(h^0 -> c c_bar)
     8.99272235E-02    2          21        21   # BR(h^0 -> g g)
     3.40982295E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.46232347E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.46232347E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.46232347E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     5.11813215E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     5.11813215E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.46232347E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.46232347E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.46232347E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     5.11813215E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     5.11813215E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     1.24226855E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     1.24226855E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     1.06480161E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     7.09867740E-03    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     4.96907418E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     4.96907418E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     4.96907418E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     4.96907418E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     4.96907418E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35     7.42785103E+00   # H^0
#    BR                NDA      ID1      ID2
     5.30468464E-04    2          13       -13   # BR(H^0 -> mu^- mu^+)
     1.49891809E-01    2          15       -15   # BR(H^0 -> tau^- tau^+)
     3.43951834E-04    2           3        -3   # BR(H^0 -> s s_bar)
     8.00273900E-01    2           5        -5   # BR(H^0 -> b b_bar)
     4.61066759E-02    2           6        -6   # BR(H^0 -> t t_bar)
     8.84050801E-04    2          23        23   # BR(H^0 -> Z Z)
     1.67844263E-03    2          24       -24   # BR(H^0 -> W^+ W^-)
     1.01790344E-04    2          25        25   # BR(H^0 -> h^0 h^0)
     1.86254782E-04    2          21        21   # BR(H^0 -> g g)
DECAY        36     7.41090815E+00   # A^0
#    BR                NDA      ID1      ID2
     5.31830228E-04    2          13       -13   # BR(A^0 -> mu^- mu^+)
     1.50277111E-01    2          15       -15   # BR(A^0 -> tau^- tau^+)
     3.44835590E-04    2           3        -3   # BR(A^0 -> s s_bar)
     8.02330144E-01    2           5        -5   # BR(A^0 -> b b_bar)
     4.64200203E-02    2           6        -6   # BR(A^0 -> t t_bar)
DECAY        37     6.55575220E+00   # H^+
#    BR                NDA      ID1      ID2
     6.01720084E-04    2         -13        12   # BR(H^+ -> mu^+ nu_e)
     1.70025604E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     3.43178580E-04    2          -3         4   # BR(H^+ -> s_bar c)
     8.28947723E-01    2          -5         6   # BR(H^+ -> b_bar t)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.01630357E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    2.65673696E+02    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    2.65690000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99938637E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    3.82514801E-03    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    3.76378486E-03        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.01630357E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    2.65673696E+02    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    2.65690000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    9.99999752E-01        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    2.48098910E-07        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    9.99999752E-01        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    2.48098910E-07        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15133433E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    2.81566322E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.04067577E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    2.48098910E-07        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    9.99999752E-01        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.02076254E-04   # BR(b -> s gamma)
    2    1.59199908E-06   # BR(b -> s mu+ mu-)
    3    3.54409066E-05   # BR(b -> s nu nu)
    4    5.19086566E-11   # BR(Bd -> mu+ mu-)
    5    3.39124654E-09   # BR(Bs -> mu+ mu-)
    6    1.10288432E-04   # BR(B_u -> tau nu)
    7    9.97277324E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73956162E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.03308169E+01   # |Delta(M_Bs)| [ps^-1] 
   20    1.42003664E-14   # Delta(g-2)_electron/2
   21    4.06996403E-11   # Delta(g-2)_muon/2
   22   -6.13753698E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    3.62695232E-03   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87801619E-01   # C7
     0305 4422   00   1    -3.22793085E-03   # C7
     0305 4322   00   1    -6.98981889E-05   # C7'
     0305 6421   00   0    -9.49035058E-02   # C8
     0305 6421   00   1    -7.82622491E-03   # C8
     0305 6321   00   1    -1.65226424E-04   # C8'
 03051111 4133   00   0    -4.84404871E-01   # C9 e+e-
 03051111 4133   00   1     1.61988002E-04   # C9 e+e-
 03051111 4233   00   1    -3.55093715E-06   # C9' e+e-
 03051111 4137   00   0    -3.92558420E+00   # C10 e+e-
 03051111 4137   00   1     1.03378831E-03   # C10 e+e-
 03051111 4237   00   1     1.08157907E-04   # C10' e+e-
 03051313 4133   00   0    -4.84404871E-01   # C9 mu+mu-
 03051313 4133   00   1     1.61977832E-04   # C9 mu+mu-
 03051313 4233   00   1    -3.55098311E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.92558420E+00   # C10 mu+mu-
 03051313 4137   00   1     1.03377813E-03   # C10 mu+mu-
 03051313 4237   00   1     1.08157953E-04   # C10' mu+mu-
 03051212 4237   00   0     1.49077549E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -2.50792468E-04   # C11 nu_1 nu_1
 03051212 4137   00   1    -2.61555702E-05   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49077549E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -2.50792468E-04   # C11 nu_2 nu_2
 03051414 4137   00   1    -2.61555594E-05   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49077549E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -2.50792459E-04   # C11 nu_3 nu_3
 03051616 4137   00   1    -2.61525141E-05   # C11' nu_3 nu_3
