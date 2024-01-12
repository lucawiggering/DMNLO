* model.h

#include "../util/DMNLO_Model.h"

* here one can include process-specific parameters

	double precision f_T(2,2,3), Delta_q(2,2,3)
	double precision m_Nucleon(2), m_Proton, m_Neutron, Mf_high(4,3), Mf_low(4,3), mu_Nucleon(2)
	complex*16 g_4Fa(2), g_4Fs(2), g_4Fap, g_4Fan, g_4Fsp, g_4Fsn
	complex*16 g_4FaContribs(2,9), g_4FsContribs(2,9)
	double precision sigmaXa(2), sigmaXs(2), sigmaXPa, sigmaXPs, sigmaXNa, sigmaXNs
	double precision sigmaXaContribs(2,9), sigmaXsContribs(2,9)
	double precision s_Mandel(4,3), t_Mandel(4,3), u_Mandel(4,3)
	integer Supp_Mode, Loop_Mode, Denom_Mode, Box_Mode, Nucl_Mode
	double precision pbp2(4,3), pap2(4,3), p1p2(4,3), papb(4,3)
	double precision Running_SI(4,3), Running_SD(4,3), Q_low

	common /Switches/ Supp_Mode, Loop_Mode, Denom_Mode, Box_Mode, Nucl_Mode
	common /Nucl_factors/ f_T, Delta_q
	common /Nucl_Masses/ m_Nucleon, m_Proton, m_Neutron, Mf_high, Mf_low, mu_Nucleon
	common /FourFermionCouplings/ g_4Fa, g_4Fs, g_4Fap, g_4Fan, g_4Fsp, g_4Fsn
	common /FourFermionCouplings/ g_4FaContribs, g_4FsContribs
	common /Sigmas/ sigmaXa, sigmaXs, sigmaXPa, sigmaXPs, sigmaXNa, sigmaXNs
	common /Sigmas/ sigmaXaContribs, sigmaXsContribs
	common /Mandelstam/ s_Mandel, t_Mandel, u_Mandel
	common /FourMomenta/ pbp2, pap2, p1p2, papb
	common /Running/Running_SI,Running_SD, Q_low