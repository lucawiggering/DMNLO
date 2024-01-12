import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm

pcmMin = 8
pcmMax = 1200
N = 25
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO = np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --renscheme 1 --particleA 1000022 --particleB 1000006 --particle1 6 --particle2 21 --pcm {} {}'

for i in tqdm(range(N)):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    print(out)
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

np.savetxt('test.txt', np.swapaxes(np.array([pcm,vsLO,vsNLO,vsNLO-vsLO]),0,1))

plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='blue',label=r'$\sigma v^{\rm NLO}$')

pcm,tree,calchep,dNLO,error = np.loadtxt('reference/Fig7_vs.txt',unpack=True)

plt.plot(pcm,tree,label='LO check')
plt.plot(pcm,tree+dNLO,label='NLO check')

plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.legend()
plt.savefig('Fig7.pdf')