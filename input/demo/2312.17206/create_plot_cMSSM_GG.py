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
N = 150
file = 'DMNLO.in'

pcm = np.linspace(pcmMin,pcmMax,N)

muR = 1368.2
command = './../../../dmnlo {} --renscheme 1 --slha {} --muR {} --particleA 1000006 --particleB -1000006 --particle1 21 --particle2 21 --pcm {} {}'


mneut = 606.295
mstop = 648.338 

def calculate_sigmav_tree(pcm,alpha_s):
    result = 0.0
    s = 4 * (mstop**2 + pcm**2)

    rho = 4 * mstop**2 / s
    beta = np.sqrt(1 - rho)

    result = (16.0 * np.pi * alpha_s**2) / (27.0 * s * beta) * (beta * (1 + rho) + rho * (rho - 2) * np.arctanh(beta))

    result += (5.0 / 2.0) * (16.0 * np.pi * alpha_s**2) / (27.0 * s * beta) * (beta * (1 + rho) + rho * (rho - 2) * np.arctanh(beta))

    result += (8.0 * np.pi * alpha_s**2) / (9.0 * s * beta) * (beta * (1 + 8 * rho) - 3 * rho * (rho + 2) * np.arctanh(beta))

#    v = 2 * beta
#    result /= v

    return result

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
        print(out)
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
ax2.get_yaxis().set_visible(False)
ax1.get_xaxis().set_visible(False)

ax1.set_zorder(ax2.get_zorder()+1)
ax1.patch.set_visible(False)

ax3.set_zorder(ax1.get_zorder()+1)
ax3.patch.set_visible(False)


T = mneut/27.5
sqrtS = np.sqrt(sman(pcm,mstop,mstop))
sigma= sqrtS*(pcm/sqrtS)**2*kv(1,sqrtS/T)*3e18
print(sigma)

ax2.plot(pcm,sigma,'-',linewidth=1.2,color='grey')
ax2.set_ylim(min(sigma),max(sigma)*1.1)
ax2.fill_between(pcm,0,sigma, facecolor='lightgrey', interpolate=True, alpha=0.8)

#ax1.set_ylim(min(vsLO*0.4),6e-7)
ax1.set_ylim(min(vsLO*0.4)*1e8,9)

ax1.fill_between(pcm,vsLO_half_mu*1e8,vsLO_2mu*1e8,color='lightgreen',alpha=0.8)
ax1.fill_between(pcm,vsNLO_half_mu*1e8,vsNLO_2mu*1e8,color='lightblue',alpha=0.8)
ax1.fill_between(pcm,vsFull_half_mu*1e8,vsFull_2mu*1e8,color='lightcoral',alpha=0.8)

ax1.plot(pcm,vsLO*1e8,color='black',ls='dashed',label=r'$\sigma^{\rm Tree}$ $(\sigma^{\rm Tree}/\sigma^{\rm MO})$')
ax1.plot(pcm,vsNLO*1e8,color=myblue,label=r'$\sigma^{\rm NLO}$ $(\sigma^{\rm NLO}/\sigma^{\rm Tree})$')
ax1.plot(pcm,vsFull*1e8,color=myred,label=r'$\sigma^{\rm Full}$ $(\sigma^{\rm Full}/\sigma^{\rm Tree})$')

# alphaS =  9.776791e-02
#ax1.plot(pcm,calculate_sigmav_tree(pcm,9.776791e-02),label='TEST',ls=':')
ax1.set_xlim(pcm[0],pcm[-1])

pcm_MO,s_MO_GG = np.loadtxt('cMSSM_MO_GG.txt',unpack=True)
pcm_MO,s_MO_uU = np.loadtxt('cMSSM_MO_uU.txt',unpack=True)
vs_MO = (s_MO_GG + 4*s_MO_uU)*vrel(pcm_MO,mstop,mstop)
ax1.plot(pcm_MO,vs_MO*1e8,label=r'$\sigma^{\mathrm{MO}}$ $(\sigma^{\rm NLO}/\sigma^{\rm MO})$',color=myorange)

plt.figtext(0.718, 0.3,r'$\tilde{t}_1 \tilde{t}^\ast_1\to g g, q \bar{q}$', fontsize=15, color=('black'))
ax3.set_xlabel(r'$p_{\rm cm}$ $\mathrm{[GeV]}$',fontsize=14)
ax1.set_ylabel(r'$\sigma v$ $[10^{-8} \, \mathrm{GeV}^{-2}]$',fontsize=14)
ax1.legend(loc='upper right',frameon=False,fontsize=13,ncol=1)

ax3.plot(pcm,vsNLO/vsLO,color=myblue)
ax3.plot(pcm,vsLO/np.interp(pcm,pcm_MO,vs_MO),color='black',ls='--')
ax3.plot(pcm,vsNLO/np.interp(pcm,pcm_MO,vs_MO),color=myorange)
ax3.plot(pcm,vsFull/vsLO,color=myred)
ax3.set_ylim(0.4,2)
#ax3.set_yticks([1,3,5])
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


plt.savefig('cMSSM_GG.pdf')