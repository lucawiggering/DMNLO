import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm

pcmMin = 4
pcmMax = 900
N = 100
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)
vsLO,vsNLO,vsSomm,vsFull = np.empty(N),np.empty(N),np.empty(N),np.empty(N)
command  = './../../../dmnlo {} --renscheme 1 --slha Scenario.spc --particleA 1000006 --particleB -1000006 --particle1 21 --particle2 21 --pcm {} {}'

for i in tqdm(range(N)):
    out = sp.getoutput(command.format('--lo  ',pcm[i],file))
    out = out.splitlines()[-1].split()[-1]
    vsLO[i] = float(out)

    out = sp.getoutput(command.format('--nlo',pcm[i],file)).splitlines()[-1].split()[-1]
    vsNLO[i] = float(out)

    out = sp.getoutput(command.format('--sommerfeld',pcm[i],file)).splitlines()[-1].split()[-1]
    print(command.format('--sommerfeld',pcm[i],file))
    vsSomm[i] = float(out)

    out = sp.getoutput(command.format('--full',pcm[i],file)).splitlines()[-1].split()[-1]
    vsFull[i] = float(out)


plt.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma v^{\rm Tree}$')
plt.plot(pcm,vsNLO,color='orange',label=r'$\sigma v^{\rm NLO}$')
plt.plot(pcm,vsFull,color='blue',label=r'$\sigma v^{\rm Full}$')
plt.plot(pcm,vsSomm,color='dodgerblue',label=r'$\sigma v^{\rm Som}$')


pcm, calchepgg, tree_ggS, tree_gg8S, tree_gg8A,calchepqqbar,treeqqbar,correction_ghost,resum,NLO_subtraction=np.loadtxt('reference/stsT2xx_Xsec_ghosts.txt',unpack=True)
pcm,sommerfeld = np.loadtxt('reference/stsT2xx_Xsec_Sommerfeld.txt',unpack=True)

tree = tree_ggS+tree_gg8S+tree_gg8A+treeqqbar
plt.plot(pcm,tree,linewidth=1.2,color='purple',alpha=0.3,label=r'$(\sigma v)^{\rm Tree}$ pub')
plt.plot(pcm,tree+correction_ghost,ls='--',alpha=0.7,linewidth=1.2,color='green',label=r'$(\sigma v)^{\rm NLO}$ pub')
plt.plot(pcm,sommerfeld,'--',linewidth=1.2,color='blue',label=r'$(\sigma v)^{\rm Som}$ pub')
plt.plot(pcm,tree+correction_ghost+resum,ls='--',alpha=0.8,linewidth=1.2,color='red',label=r'$(\sigma v)^{\rm Full}_{}$ pub')

plt.xlabel(r'$p_{\rm cm}$ [GeV]')
plt.legend()
plt.yscale('log')
plt.savefig('Fig16a.pdf')