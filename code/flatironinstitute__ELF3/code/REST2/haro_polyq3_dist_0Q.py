import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_context('talk')

q290=np.loadtxt("haro_pq3_290K.xvg")
q300=np.loadtxt("haro_pq3_300K.xvg")
q320=np.loadtxt("haro_pq3_320K.xvg")
q405=np.loadtxt("haro_pq3_405K.xvg")

N_points=50001
n_bins=250

plt.figure(figsize=(8.5,6), dpi=300)
plt.hist(q290,n_bins,label=['290K'], color='blue')
plt.hist(q300,n_bins,label=['300K'], color='green',alpha=0.8)
plt.hist(q320,n_bins,label=['320K'], color='orange',alpha=0.8)
plt.hist(q405,n_bins,label=['405K'], color='red',alpha=0.8)
plt.legend(loc='upper right')
plt.xlabel(r'H$\rm_{aro}$-PolyQ3 Distance (nm)', fontsize=22)
plt.ylabel('# Structures', fontsize=22)
#plt.text(s="320K",x=1.5,y=8000,fontsize=22)
plt.text(s="0Q",x=4,y=6100,fontsize=24)
plt.xlim([0,8])
#plt.ylim([0,10000])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('fig2.png',bbox_inches='tight')
