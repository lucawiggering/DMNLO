import pyslha
import subprocess as sp
import numpy as np
import parse 
import matplotlib.pyplot as plt

# path to the SPheno binary
SPHENO_PATH = '../../../external/SPheno-3.2.3/bin/SPheno'

# SPheno input file
DEFAULT_INPUT = 'ScenarioA.in'
doc = pyslha.read(DEFAULT_INPUT,ignorenomass=True) # save the default input file in doc

# new SPheno input file
SPHENO_INPUT_FILE = 'SPheno.in'

# SPheno output file
SPHENO_OUTPUT_FILE = 'LesHout.spc'

# spcification of the region of the parameter scan
N = 100
M1 = np.linspace(100,1000,N)
SIpLO,SInLO,SDpLO,SDnLO = np.empty(N),np.empty(N),np.empty(N),np.empty(N)
SIpNLO,SInNLO,SDpNLO,SDnNLO = np.empty(N),np.empty(N),np.empty(N),np.empty(N)

# command to run DM@NLO in direct detection mode
command  = './../../../dmnlo {} --DD --slha LesHout.spc --legacy --renscheme 1 DMNLO.in'


for i in range(N): # loop over M1

    # modify soft-breaking parameters 
    doc.blocks['EXTPAR'][1] = M1[i] # M1

    # write new SPheno input
    pyslha.write(SPHENO_INPUT_FILE,doc)

    # run SPheno to compute the spectrum
    sp.check_call([SPHENO_PATH,SPHENO_INPUT_FILE, SPHENO_OUTPUT_FILE])

    # leading order results
    out = sp.getoutput(command.format('--lo')).splitlines()
    SI_line, SD_line = out[-3],out[-1]
    SIpLO[i], SInLO[i] = parse.parse('proton SI : {:e}   neutron SI : {:e}',SI_line)
    SDpLO[i], SDnLO[i] = parse.parse('proton SD : {:e}   neutron SD : {:e}',SD_line)

    # next-to-leading order results
    out = sp.getoutput(command.format('--nlo')).splitlines()
    SI_line, SD_line = out[-3],out[-1]
    SIpNLO[i], SInNLO[i] = parse.parse('proton SI : {:e}   neutron SI : {:e}',SI_line)
    SDpNLO[i], SDnNLO[i] = parse.parse('proton SD : {:e}   neutron SD : {:e}',SD_line)


# create plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
fig.tight_layout(pad=3.0)




ax1.plot(M1,SIpLO,color='black',label='Tree')
ax1.plot(M1,SIpNLO,color='royalblue',label='NLO')
ax1.set_xlabel(r'$M_1$ [GeV]')
ax1.set_ylabel(r'$\sigma_p^{\rm SI}$ [$\mathrm{cm}^2$]')
ax1.set_yscale('log')
ax1.legend(frameon=False)

M1ref,SIpLOref = np.loadtxt('reference/Fig5a_LO.txt',unpack=True,delimiter=',')
ax1.plot(M1ref,SIpLOref,color='red',ls='dotted')
M1ref,SIpNLOref = np.loadtxt('reference/Fig5a_NLO.txt',unpack=True,delimiter=',')
ax1.plot(M1ref,SIpNLOref,color='orange',ls='dotted')

ax2.plot(M1,SInLO,color='black',label='Tree')
ax2.plot(M1,SInNLO,color='royalblue',label='NLO')
ax2.set_xlabel(r'$M_1$ [GeV]')
ax2.set_ylabel(r'$\sigma_n^{\rm SI}$ [$\mathrm{cm}^2$]')
ax2.set_yscale('log')
ax2.legend(frameon=False)


M1ref,SIpLOref = np.loadtxt('reference/Fig5b_LO.txt',unpack=True,delimiter=',')
ax2.plot(M1ref,SIpLOref,color='red',ls='dotted')
M1ref,SIpNLOref = np.loadtxt('reference/Fig5b_NLO.txt',unpack=True,delimiter=',')
ax2.plot(M1ref,SIpNLOref,color='orange',ls='dotted')

ax3.plot(M1,SDpLO,color='black',label='Tree')
ax3.plot(M1,SDpNLO,color='royalblue',label='NLO')
ax3.set_xlabel(r'$M_1$ [GeV]')
ax3.set_ylabel(r'$\sigma_p^{\rm SD}$ [$\mathrm{cm}^2$]')
ax3.set_yscale('log')
ax3.legend(frameon=False)

M1ref,SIpNLOref = np.loadtxt('reference/Fig5c_NLO.txt',unpack=True,delimiter=',')
ax3.plot(M1ref,SIpNLOref,color='orange',ls='dotted')

ax4.plot(M1,SDnLO,color='black',label='Tree')
ax4.plot(M1,SDnNLO,color='royalblue',label='NLO')
ax4.set_xlabel(r'$M_1$ [GeV]')
ax4.set_ylabel(r'$\sigma_n^{\rm SD}$ [$\mathrm{cm}^2$]')
ax4.set_yscale('log')
ax4.legend(frameon=False)

M1ref,SIpNLOref = np.loadtxt('reference/Fig5d_NLO.txt',unpack=True,delimiter=',')
ax4.plot(M1ref,SIpNLOref,color='orange',ls='dotted')


plt.savefig('Fig5_legacy.pdf')