import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("talk")

q0=np.loadtxt('0Q_290K.pc1')
q7=np.loadtxt('7Q_290K.pc1')
q13=np.loadtxt('13Q_290K.pc1')
q19=np.loadtxt('19Q_290K.pc1')

hq0=np.loadtxt('hightemp0q.dat')
hq7=np.loadtxt('higihtemp7q.dat')
hq13=np.loadtxt('hightemp13q.dat')
hq19=np.loadtxt('hightemp19q.dat')

fig, axs = plt.subplots(1,2,figsize=(7,1.75), dpi=300)
#plt.subplots_adjust(hspace=0)
axs[0].plot(-q7, color='green')
axs[0].plot(q13, color='gold')
axs[0].plot(q19, color='deepskyblue')
axs[0].plot(q0, color='magenta')

axs[1].plot(hq7,color='cyan')
axs[1].plot(hq13,color='darkblue')
axs[1].plot(-hq19,color='skyblue')
axs[1].yaxis.tick_right()
axs[1].plot(hq0,color='magenta')

axs[0].set_ylim([-.3,.75])
axs[1].set_ylim([-.3,.75])

axs[0].text(s='290K',x=1700,y=.50,fontsize=14, weight='bold')
axs[1].text(s='415K',x=1875,y=.50,fontsize=14, weight='bold')

leg=axs[0].legend(['7Q', '13Q', '19Q', '0Q'],fontsize=10, loc="upper left")
x=range(186)
axs[0].tick_params(axis='x')
axs[0].tick_params(axis='y')
axs[1].tick_params(axis='x')
axs[1].tick_params(axis='y')
leg.get_frame().set_alpha(0)

#axs[0,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[1,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[2,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[3,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[3,0].set_xticklabels([0,500,1000,1500,2000,2500])

axs[0].set_xlim(1250,1850)
axs[1].set_xlim(1300,2100)

#axs[3,1].set_xticks([1400,1800,2200])
#axs[3,1].set_xticklabels(['1400','1800','2200'])

axs[0].set_xticks([1400,1700])
axs[0].set_xticklabels(["100","120"])
axs[1].set_xticks([1500,1900])
axs[1].set_xticklabels(["90","110"])
axs[0].set_yticks([-.3,0,0.3,0.6])
axs[1].set_yticks([-.3,0,0.3,0.6])

axs[0].set_ylabel("PC1",labelpad=.5)
axs[0].set_xlabel("Approx. Residue Index")
axs[0].xaxis.set_label_coords(.95, -.3)
axs[0].yaxis.set_label_coords(-.2,.5)
#axs[0,0].get_xaxis().set_visible(False)
#axs[1,0].get_xaxis().set_visible(False)
#axs[2,0].get_xaxis().set_visible(False)

#axs[3,0].text(2000, -.6, 'Contact Index', ha='center', va='center',fontsize=22)
#axs[2,0].text(775,.7, 'PC1', ha='center', va='center',fontsize=22,rotation=90)
#axs[3,0].yaxis.set_label_coords(-.1,2.25)
plt.subplots_adjust(hspace=0,wspace=0)

#axs[1].axvspan(112,118,color='deepskyblue',alpha=0.5, lw=0)
#axs[2].axvspan(112,124,color='deepskyblue',alpha=0.5, lw=0)
#axs[3].axvspan(112,131,color='deepskyblue',alpha=0.5, lw=0)


plt.savefig('ipca_2col_290K_415K.png',bbox_inches='tight')
