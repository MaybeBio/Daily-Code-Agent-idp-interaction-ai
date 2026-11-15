import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("talk")

q0=np.loadtxt('0q_290K.pc1')
q7=np.loadtxt('7q_290K.pc1')
q13=np.loadtxt('13q_290K.pc1')
q19=np.loadtxt('19q_290K.pc1')

hq0=np.loadtxt('0q_290K.pc2')
hq7=np.loadtxt('7q_290K.pc2')
hq13=np.loadtxt('13q_290K.pc2')
hq19=np.loadtxt('19q_290K.pc2')

fig, axs = plt.subplots(4,2,figsize=(10,8), dpi=300)
#plt.subplots_adjust(hspace=0)
axs[0,0].plot(q0, color='blue')
axs[0,1].plot(hq0,color='red')
axs[0,1].yaxis.tick_right()

axs[1,0].plot(-q7, color='blue')
axs[1,1].plot(hq7,color='red')
axs[1,1].yaxis.tick_right()
axs[2,1].yaxis.tick_right()
axs[3,1].yaxis.tick_right()

axs[2,0].plot(q13, color='blue')
axs[2,1].plot(hq13,color='red')

axs[3,0].plot(q19, color='blue')
axs[3,1].plot(-hq19,color='red')

axs[0,0].set_ylim([-.6,.6])
axs[0,1].set_ylim([-.6,.6])
axs[1,0].set_ylim([-.6,.6])
axs[1,1].set_ylim([-.6,.6])
axs[2,0].set_ylim([-.6,.6])
axs[2,1].set_ylim([-.6,.6])
axs[3,0].set_ylim([-.6,.6])
axs[3,1].set_ylim([-.6,.6])

axs[0,0].text(s='0Q',x=200,y=.35, fontsize=22)
axs[0,0].text(s='320K',x=1750,y=.35,fontsize=22)
axs[1,0].text(s='7Q',x=200,y=.31,fontsize=22)
axs[2,0].text(s='13Q',x=200,y=.29,fontsize=22)
axs[3,0].text(s='19Q',x=200,y=.28,fontsize=22)


x=range(186)
axs[0,0].tick_params(axis='x')
axs[0,0].tick_params(axis='y')
axs[1,0].tick_params(axis='x')
axs[1,0].tick_params(axis='y')
axs[2,0].tick_params(axis='x')
axs[2,0].tick_params(axis='y')
axs[3,0].tick_params(axis='x')
axs[3,0].tick_params(axis='y')

#axs[0,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[1,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[2,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[3,0].set_xticks([0,500,1000,1500,2000,2500])
#axs[3,0].set_xticklabels([0,500,1000,1500,2000,2500])

axs[0,0].set_xlim(0,2320)
axs[1,0].set_xlim(0,2340)
axs[2,0].set_xlim(0,2480)
axs[3,0].set_xlim(0,2500)
axs[0,1].set_xlim(0,2360)
axs[1,1].set_xlim(0,2410)
axs[2,1].set_xlim(0,2450)
axs[3,1].set_xlim(0,2500)

#axs[3,0].set_xticks([1100,1500,1900])
#axs[3,0].set_xticklabels(['1100','1500','1900'])
axs[0,1].set_xticks([])
axs[1,1].set_xticks([])
axs[2,1].set_xticks([])
axs[0,0].set_xticks([])
axs[1,0].set_xticks([])
axs[2,0].set_xticks([])
axs[3,0].set_xticks([0,1000,2000])
axs[3,0].set_xticklabels(['0','1000','2000'])
#axs[3,0].set_xticks([1000,1500,2000],['1000','1500','2000'])

axs[0,0].set_yticks([-0.6,-0.3,0,0.3,0.6])
axs[1,0].set_yticks([-0.6,-0.3,0,0.3,0.6])
axs[2,0].set_yticks([-0.6,-0.3,0,0.3,0.6])
axs[3,0].set_yticks([-0.6,-0.3,0,0.3,0.6])

axs[0,1].set_yticks([-0.3,0,0.3,0.6])
axs[1,1].set_yticks([-0.3,0,0.3,0.6])
axs[2,1].set_yticks([-0.3,0,0.3,0.6])
axs[3,1].set_yticks([-0.3,0,0.3,0.6])

axs[0,0].get_xaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[2,0].get_xaxis().set_visible(False)


#axs[3,0].set_xlabel('Contact Index',fontsize=24)
#axs[3,0].xaxis.set_label_coords(2000,-.3)
axs[3,0].text(2500, -1.05, 'Contact Index', ha='center', va='center',fontsize=22)
#axs[2,0].set_ylabel('PC1',fontsize=24)
axs[2,0].text(-600,.7, 'PC1', ha='center', va='center',fontsize=22,rotation=90)
axs[2,1].text(3000,.7, 'PC2', ha='center', va='center',fontsize=22,rotation=90)
axs[3,0].yaxis.set_label_coords(-.1,2.25)
plt.subplots_adjust(hspace=0,wspace=0)

#axs[1].axvspan(112,118,color='deepskyblue',alpha=0.5, lw=0)
#axs[2].axvspan(112,124,color='deepskyblue',alpha=0.5, lw=0)
#axs[3].axvspan(112,131,color='deepskyblue',alpha=0.5, lw=0)


plt.savefig('pc1_2_290K.png',bbox_inches='tight')
