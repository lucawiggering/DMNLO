import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm
import concurrent.futures 

pcmMin = 5
pcmMax = 1200
N = 45
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO,vsSomm,vsFull = np.empty(N),np.empty(N),np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --slha ScenarioI.spc --renscheme 1 --muR 1784.64 --particleA 1000006 --particleB 1000006 --particle1 6 --particle2 6 --pcm {} {}'

def compute(i):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

    out = sp.getoutput(command.format('--sommerfeld',pcm[i],file)).splitlines()[-1].split()[-1]
    vsSomm[i] = float(out)

    out = sp.getoutput(command.format('',pcm[i],file)).splitlines()[-1].split()[-1]
    vsFull[i] = float(out)

pcmDip,tree,calchep,dNLO_Dipole,errorDip= np.loadtxt('reference/Dipole_tt_wError.txt',unpack=True)
plt.plot(pcmDip,tree+dNLO_Dipole,color='black',ls='dotted',alpha=0.8,label=r'$\sigma^{\rm NLO} dipole$')

with concurrent.futures.ThreadPoolExecutor(max_workers = 7) as executor:
    results = [executor.submit(compute,i) for i in range(N)]

plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v_{\rm Tree}$')
plt.plot(pcm,vsNLO,color='green',label=r'$\sigma v^{\rm NLO}$')
plt.plot(pcm,vsSomm,color='purple',ls='dashdot',label=r'$\sigma v_{\rm Som}$')
plt.plot(pcm,vsFull,color='blue',label=r'$\sigma v_{\rm F}$')

pcmCheck, vsLOCheck = np.loadtxt('reference/Fig10a_vsLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck,color='orange',label='LO pub',alpha=0.5)
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig10a_vsNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck,color='gray',label='NLO pub',alpha=0.8,ls='dashed')
pcmCheck, vsFullCheck = np.loadtxt('reference/Fig10a_vsFull.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsFullCheck,color='hotpink',label='Full pub',alpha=0.8,ls='dashed')




plt.title(r'$\tilde{t}_1 \tilde{t}_1 \to t t$')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.ylim(8e-10,9e-9)
plt.yscale('log')
plt.legend(frameon=False)
plt.savefig('Fig11a.pdf')