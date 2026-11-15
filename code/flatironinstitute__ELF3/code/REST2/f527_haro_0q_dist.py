import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_context('talk')

q7=np.loadtxt("7q_haro_f96_290K.dat")
q0=np.loadtxt("0q_haro_f96_290K.dat")
N_points=50001
n_bins=250

plt.figure(figsize=(8.5,6), dpi=300)
plt.hist([q7],n_bins,label=['7Q'], color='C0')
plt.hist([q0],n_bins,label=['0Q'], color='C2',alpha=0.8)
plt.legend(loc='upper right')
plt.xlabel(r'F527-H$_{\rm aro}$ Distance (nm)', fontsize=22)
plt.ylabel('# Structures', fontsize=22)
plt.text(s="290K",x=1.5,y=8000,fontsize=22)
plt.xlim([0,3])
plt.ylim([0,10000])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('haro_f96_0Q.png',bbox_inches='tight')
