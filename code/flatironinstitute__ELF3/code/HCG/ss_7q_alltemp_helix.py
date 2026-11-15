import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

q0_290=np.loadtxt('0q_HG290_norm')
q0_300=np.loadtxt('0q_HG300_norm')
q0_320=np.loadtxt('0q_HG320_norm')
q0_415=np.loadtxt('0q_HG415_norm')

q7_290=np.loadtxt('7q_HG290_norm')
q7_300=np.loadtxt('7q_HG300_norm')
q7_320=np.loadtxt('7q_HG320_norm')
q7_415=np.loadtxt('7q_HG415_norm')

default_ticks=range(432,605)

#plt.plot(q0_290)
#plt.plot(q0_300)
fig = plt.figure(figsize=(9,3), dpi=300)
#fig = plt.figure(figsize=(8.5,6), dpi=300)
plt.plot(default_ticks,q7_290, color='blue')
plt.plot(default_ticks,q7_300, color='green')
plt.plot(default_ticks,q7_320, color='orange')
plt.plot(default_ticks,q7_415, color='red')
#plt.xlim([0,173])
plt.ylim([0,0.34])


plt.xticks(default_ticks[::20],default_ticks[::20])
plt.tick_params(axis='x',labelsize=18)
plt.tick_params(axis='y',labelsize=18)
plt.yticks([0.0, 0.1, 0.2, 0.3],["0","10", "20", "30"])
#plt.axvspan(494,504,color='yellow',alpha=0.5, lw=0)
plt.axvspan(544,550,color='deepskyblue',alpha=0.5, lw=0)
plt.axvline(x=527,color='red',ymin=0,ymax=1,ls='--',alpha=0.5)
#plt.xticklabels(default_ticks[::20])
plt.legend(['290K','300K','320K','415K'],loc='upper left',fontsize=14)
plt.xlabel('Residue Index',fontsize=24)
plt.ylabel('Percent %',fontsize=24)

#plt.subplots_adjust(hspace=0)

plt.savefig('helicity_per_residue_hcg.png',bbox_inches='tight')
