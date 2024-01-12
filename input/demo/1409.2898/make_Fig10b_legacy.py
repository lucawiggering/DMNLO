import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm

pcmMin = 5
pcmMax = 500
N = 40
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO = np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --renscheme 1 --legacy --slha Scenario.spc --muR 1000 --particleA 1000022 --particleB 1000006 --particle1 6 --particle2 25 --pcm {} {}'


for i in tqdm(range(N)):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

pcmCheck, vsLOCheck = np.loadtxt('reference/Fig10b_vsLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck*1e-8,color='green',label='LO pub')
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig10b_vsNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck*1e-8,color='purple',label='NLO pub',ls='dashed')


plt.title(r'$\tilde{\chi}\tilde{t}\to t h$ legacy')
plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='blue',label=r'$\sigma v^{\rm NLO}$')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.legend(frameon=False)
plt.savefig('Fig10b_legacy.pdf')