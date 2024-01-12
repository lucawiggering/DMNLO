import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm

pcmMin = 5
pcmMax = 500
N = 25
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO = np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --renscheme 1 --slha Scenario.spc --muR 1000 --particleA 1000022 --particleB 1000006 --particle1 5 --particle2 24 --pcm {} {}'


for i in tqdm(range(N)):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    print(out)
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

pcmCheck, vsLOCheck = np.loadtxt('reference/Fig10c_vsLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck,color='green',label='LO pub')
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig10c_vsNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck,color='purple',label='NLO pub')

for i in range(N):
    print(pcm[i],vsNLO[i])

plt.title(r'$\tilde{\chi}\tilde{t}\to b W^+$')
plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='blue',label=r'$\sigma v^{\rm NLO}$')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.legend(frameon=False)
plt.savefig('Fig10c.pdf')