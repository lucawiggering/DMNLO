import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from tqdm import tqdm
import matplotlib.gridspec as gridspec
from scipy.special import kv


plt.rc('text', usetex=True)

myred='#DA0002'
myblue ='#006DDC'
myorange ='#FF8000'
mygreen ='#01A801'

pcmMin = 5
pcmMax = 400
N = 10
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)

muR = 1368.2
command = './../../../dmnlo {} --renscheme 1 --slha {} --muR {} --particleA 1000006 --particleB -1000006 --particle1 22 --particle2 21 --pcm {} {}'


mneut = 606.295
mstop = 648.338 

def Kallen(x,y,z):
    return x**2 + y**2 + z**2 - 2*x*y - 2*y*z - 2*x*z

def sman(pcm,ma,mb):
    return (np.sqrt(pcm**2 + ma**2) + np.sqrt(pcm**2 + mb**2))**2

def vrel(pcm,ma,mb):
    s = sman(pcm,ma,mb)
    return 2*np.sqrt(Kallen(s,ma**2,mb**2))/s

def compute(scale,slha):
    vsLO,vsNLO,vsFull = np.empty(N),np.empty(N),np.empty(N)
    
    for i in tqdm(range(N)):
        out = sp.getoutput(command.format('--lo  ',slha,scale,pcm[i],file))
        out = out.splitlines()[-1].split()[-1]
        vsLO[i] = float(out)

        out = sp.getoutput(command.format('--nlo',slha,scale,pcm[i],file)).splitlines()[-1].split()[-1]
        vsNLO[i] = float(out)

        out = sp.getoutput(command.format('',slha,scale,pcm[i],file)).splitlines()[-1].split()[-1]
        vsFull[i] = float(out)

    return vsLO,vsNLO,vsFull

vsLO_half_mu,vsNLO_half_mu, vsFull_half_mu = compute(muR/2,'cMSSM_half_mu.spc')
vsLO,vsNLO, vsFull = compute(muR,'cMSSM.spc')
vsLO_2mu,vsNLO_2mu, vsFull_2mu = compute(2*muR,'cMSSM_2mu.spc')

fig = plt.figure(layout='constrained')
#ax1 = fig.add_subplot()
fig.subplots_adjust(hspace=0.1)
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
gs.update(hspace=0.001)
ax1=fig.add_subplot(gs[0])
ax3=fig.add_subplot(gs[1],sharex=ax1)
ax2=ax1.twinx() 
ax1.get_xaxis().set_visible(False)

T = mneut/27.5
sqrtS = np.sqrt(sman(pcm,mstop,mneut))
sigma= sqrtS*(pcm/sqrtS)**2*kv(1,sqrtS/T)*3e18
print(sigma)

#ax2.plot(pcm,sigma,'-',linewidth=1.2,color='grey',zorder=0.1)
#ax2.fill_between(pcm,0,sigma, facecolor='lightgrey', interpolate=True, alpha=0.8,zorder=0.1)
ax2.set_ylim(0,max(sigma)*1.5)
ax2.set_yscale('linear')
ax2.set_yticks([])
#ax1.set_ylim(min(vsLO*0.6),9e-7)
ax1.set_zorder(ax2.get_zorder()+1)
ax1.patch.set_visible(False)

ax3.set_zorder(ax1.get_zorder()+1)
ax3.patch.set_visible(False)

ax1.fill_between(pcm,vsLO_half_mu,vsLO_2mu,color='lightgreen',alpha=0.8)
#ax1.fill_between(pcm,vsNLO_half_mu,vsNLO_2mu,color='lightblue',alpha=0.8)
#ax1.fill_between(pcm,vsFull_half_mu,vsFull_2mu,color='lightcoral',alpha=0.8)

ax1.plot(pcm,vsLO,color='black',ls='dashed',label=r'$\sigma^{\rm Tree}$ $(\sigma^{\rm Tree}/\sigma^{\rm MO})$')
ax1.plot(pcm,vsNLO,color=myblue,label=r'$\sigma^{\rm NLO}$ $(\sigma^{\rm NLO}/\sigma^{\rm Tree})$',ls=':')
#ax1.plot(pcm,vsFull,color=myred,label=r'$\sigma^{\rm Full}$ $(\sigma^{\rm Full}/\sigma^{\rm Tree})$')
ax1.set_xlim(pcm[0],pcm[-1])

pcm_MO,s_MO = np.loadtxt('cMSSM_MO_AG.txt',unpack=True)
vs_MO = s_MO*vrel(pcm_MO,mstop,mneut)
ax1.plot(pcm_MO,vs_MO,label=r'$\sigma^{\mathrm{MO}}$ $(\sigma^{\rm NLO}/\sigma^{\rm MO})$',color=myorange)

ax1.annotate(r'$\tilde{t}_1 \tilde{t}^\ast_1\to g g$',xy=(320,4e-8),fontsize=15)
ax3.set_xlabel(r'$p_{\rm cm}$ $\mathrm{[GeV]}$',fontsize=14)
ax1.set_ylabel(r'$\sigma v$ $[\mathrm{GeV}^{-2}]$',fontsize=14)
ax1.legend(loc='upper right',frameon=False,fontsize=13,ncol=2)

ax3.plot(pcm,vsNLO/vsLO,color=myblue)
ax3.plot(pcm,vsLO/np.interp(pcm,pcm_MO,vs_MO),color='black',ls='--')
ax3.plot(pcm,vsNLO/np.interp(pcm,pcm_MO,vs_MO),color=myorange)
ax3.plot(pcm,vsFull/vsLO,color=myred)
#ax3.set_ylim(0.1,6)
ax3.set_yticks([1,3,5])
ax3.set_ylabel(r'$\mathrm{ratio}$',fontsize=14)

#ax1.set_yscale('log')
ax1.tick_params(which='both',labelsize=13)  
ax1.minorticks_on()
ax1.yaxis.set_ticks_position('both')
#ax2.tick_params(which='both',labelsize=16) 
ax3.tick_params(which='both',labelsize=13) 
ax3.minorticks_on()
ax3.grid(linestyle='--',alpha=0.4)
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')


plt.savefig('cMSSM_AG.pdf')