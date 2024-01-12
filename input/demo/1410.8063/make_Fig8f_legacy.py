import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import concurrent.futures 

pcmMin = 20
pcmMax = 1000
N = 60
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO,vsSomm,vsFull = np.empty(N),np.empty(N),np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --slha ScenarioII.spc  --legacy --renscheme 1 --particleA 1000006 --particleB -1000006 --particle1 24 --particle2 -24 --pcm {} {}'

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

with concurrent.futures.ThreadPoolExecutor(max_workers = 7) as executor:
    results = [executor.submit(compute,i) for i in range(N)]

plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v_{\rm Tree}$')
plt.plot(pcm,vsNLO,color='red',label=r'$\sigma v^{\rm NLO}$')
plt.plot(pcm,vsSomm,color='purple',ls='dashdot',label=r'$\sigma v_{\rm Som}$')
plt.plot(pcm,vsFull,color='blue',label=r'$\sigma v_{\rm F}$')

pcmCheck, vsLOCheck = np.loadtxt('reference/Fig8f_vsLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsLOCheck,color='orange',label='LO pub',alpha=0.6)
pcmCheck, vsNLOCheck = np.loadtxt('reference/Fig8f_vsNLO.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsNLOCheck,color='green',label='NLO pub',ls='dashed')
pcmCheck, vsFullCheck = np.loadtxt('reference/Fig8f_vsFull.txt',unpack=True,delimiter=',')
plt.plot(pcmCheck,vsFullCheck,color='gray',label='Full pub',ls='dashed')

plt.title(r'$\tilde{t}_1 \tilde{t}^\ast_1 \to W^+ W^-$ legacy')
plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.ylabel(r'$\sigma v$ [GeV$^{-2}$]')
plt.yscale('log')
plt.legend(frameon=False)
plt.savefig('Fig8f_legacy.pdf')