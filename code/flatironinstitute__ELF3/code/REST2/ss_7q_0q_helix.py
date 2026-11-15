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

#plt.plot(q0_290)
#plt.plot(q0_300)
fig = plt.figure(figsize=(9,3), dpi=300)
#fig = plt.figure(figsize=(8.5,6), dpi=300)
plt.plot(q7_290, color='blue')
plt.plot(q7_300, color='green')
plt.plot(q7_320, color='orange')
plt.plot(q7_415, color='red')
plt.xlim([0,173])
plt.ylim([0,0.34])

plt.text(s=r'$H_N3$',x=87,y=.08,fontsize=16,color='black')
plt.text(s=r'$H_N2$',x=92,y=.15,fontsize=16,color='black')
plt.text(s=r'$H_N1$',x=97,y=.25,fontsize=16,color='black')
plt.text(s=r'$H_C1$',x=120,y=.09,fontsize=16,color='black')

x=range(173)
plt.xticks(x[::20],x[::20])
plt.tick_params(axis='x',labelsize=18)
plt.tick_params(axis='y',labelsize=18)

plt.axvspan(63,73,color='yellow',alpha=0.5, lw=0)
plt.axvspan(113,119,color='deepskyblue',alpha=0.5, lw=0)
plt.axvspan(137,141,color='deepskyblue',alpha=0.5, lw=0)
plt.axvspan(150,154,color='deepskyblue',alpha=0.5, lw=0)
plt.axvline(x=96,color='red',ymin=0,ymax=1,ls='--',alpha=0.5)
#plt.xticklabels(x[::20])
plt.legend(['290K','300K','320K','415K'],loc='upper left',fontsize=14)
plt.xlabel('Residue Index',fontsize=24)
plt.ylabel('Ratio',fontsize=24)

#plt.subplots_adjust(hspace=0)

plt.savefig('test.png',bbox_inches='tight')
