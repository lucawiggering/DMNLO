import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm
import concurrent.futures

pcmMin = 10
pcmMax = 350
N = 50
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO,vsSomm,vsFull = np.empty(N),np.empty(N),np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --slha ScenarioI.spc --legacy --renscheme 2 --muR 1000 --particleA 1000015 --particleB -1000015 --particle1 6 --particle2 -6 --pcm {} {}'

def compute(i):
    print(i)

    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

    out = sp.getoutput(command.format('--sommerfeld',pcm[i],file)).splitlines()[-1].split()[-1]
    vsSomm[i] = float(out)

    out = sp.getoutput(command.format('',pcm[i],file)).splitlines()[-1].split()[-1]
    vsFull[i] = float(out)

with concurrent.futures.ThreadPoolExecutor(max_workers = 7) as executor:
    results = [executor.submit(compute,i) for i in range(N)]

plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v_{\rm Tree}$')
plt.plot(pcm,vsNLO,color='green',label=r'$\sigma v^{\rm NLO}$')
plt.plot(pcm,vsSomm,color='purple',ls='dashdot',label=r'$\sigma v_{\rm Som}$')
plt.plot(pcm,vsFull,color='blue',label=r'$\sigma v_{\rm F}$')

pcmCheck, vsLOCheck = np.loadtxt('reference/Fig6b_vsLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck,color='orange',label='LO pub',alpha=0.5)
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig6b_vsNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck,color='gray',label='NLO pub',alpha=0.8,ls='dashed')
pcmCheck, vsFullCheck = np.loadtxt('reference/Fig6b_vsFull.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsFullCheck,color='hotpink',label='Full pub',alpha=0.8,ls='dashed')
pcmCheck, vsFullCheck = np.loadtxt('reference/Fig6b_vsSomm.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsFullCheck,color='teal',label='Somm pub',alpha=0.5,ls='-')


plt.title(r'$\tilde{\tau}_1 \tilde{\tau}^\ast_1 \to t \bar{t}$ legacy')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.ylim(3e-9,8e-9)
plt.legend(frameon=False)
plt.savefig('Fig6b_legacy.pdf')