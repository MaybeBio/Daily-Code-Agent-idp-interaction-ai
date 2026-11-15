import matplotlib.pyplot as plt
import numpy as np

wt290 = np.loadtxt("f96_aro_290K.xvg")
wt300 = np.loadtxt("f96_aro_300K.xvg")
wt320 = np.loadtxt("f96_aro_320K.xvg")
wt405 = np.loadtxt("f96_aro_405K.xvg")
n_bins=250

plt.figure(figsize=(8.5,6), dpi=300)
plt.hist([wt290],n_bins,label=['290K'], color='blue',alpha=0.5)
plt.hist([wt300],n_bins,label=['300K'], color='green',alpha=0.5)
plt.hist([wt320],n_bins,label=['320K'], color='orange',alpha=0.5)
plt.hist([wt405],n_bins,label=['405K'], color='red',alpha=0.5)
plt.legend(loc='upper right',fontsize=20)
plt.xlabel('Minimum Distance (nm)', fontsize=20)
plt.ylabel('# Structures', fontsize=20)
plt.xlim([0,3])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('f96_haro_comp.png',bbox_inches='tight')
