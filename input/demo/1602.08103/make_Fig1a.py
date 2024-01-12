import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm

pcmMin = 5
pcmMax = 500
N = 100
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO = np.empty(N),np.empty(N)
command  = './../../../src/main {} --renscheme 1 --particleA 1000022 --particleB 1000023 --particle1 5 --particle2 -5 --pcm {} {}'


for i in tqdm(range(N)):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='blue',label=r'$\sigma v^{\rm NLO}$')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.yscale('log')
plt.legend()
plt.savefig('Fig1a.pdf')