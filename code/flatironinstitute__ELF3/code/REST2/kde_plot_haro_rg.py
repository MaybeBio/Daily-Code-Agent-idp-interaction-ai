import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

sns.set_context("talk")


data1=np.loadtxt("dist.dat")
data2=np.loadtxt("rg.dat")
print(data1)
print(data2)
matplotlib.rcParams['axes.edgecolor'] = 'C1'
#matplotlib.rcParams['ytick.color'] = 'C0'
kde=sb.kdeplot(data1,data2=data2, cmap="Reds",shade=True)
#kde.set(yticklabels=[])  
kde.set(ylabel=None)
#kde.tick_params(left=False)
fig=kde.get_figure()
#plt.tick_params(bottom=False,top=False)
#plt.xticks([])
#plt.yticks([])
plt.xlim(-.25,6)
plt.ylim(1.5,6.5)
#plt.yticks([2,3,4,5,6],['2','3','4','5','6'])
fig.savefig("wt_340K_rg_hdist.png")
