import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import concurrent.futures 

pcmMin = 5
pcmMax = 500
N = 40
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO = np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --renscheme 1 --result s --legacy --slha ScenarioI.spc --muR 1000 --particleA 1000024 --particleB 1000022 --particle1 6 --particle2 -5 --pcm {} {}'


def compute(i):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    print(command.format('--lo  ',pcm[i],file))
    print(out)
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)


with concurrent.futures.ThreadPoolExecutor(max_workers = 7) as executor:
    results = [executor.submit(compute,i) for i in range(N)]


pcmCheck, vsLOCheck = np.loadtxt('reference/Fig6a_sLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck,color='green',label='LO pub')
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig6a_sNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck,color='purple',label='NLO pub')
pcmCheck, vsMOCheck = np.loadtxt('MO/Fig6a.txt',unpack=True)
plt.plot(pcmCheck,vsMOCheck,color='pink',label='MO pub')

plt.title(r'$\tilde{\chi}^0_1\tilde{\chi}^+_1\to t \bar{b}$ legacy')
plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='blue',label=r'$\sigma v^{\rm NLO}$')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma$ [GeV$^{-2}$]')
plt.yscale('log')
plt.legend(frameon=False)
plt.savefig('Fig6a_legacy.pdf')