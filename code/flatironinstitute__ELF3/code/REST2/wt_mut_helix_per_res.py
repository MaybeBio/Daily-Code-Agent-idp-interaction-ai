import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/parts/helix_500ns_nowat.xtc'
top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/helix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp1 = md.compute_dssp(traj)

boolean1_1 = dssp1 =='H'
boolean1_2 = dssp1 =='G'
zeros_and_ones1 = 1* (boolean1_1+boolean1_2)
helix_prob1 = zeros_and_ones1.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir2/parts/helix_500ns_nowat2.xtc'
#top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp2 = md.compute_dssp(traj)

boolean2_1 = dssp2 =='H'
boolean2_2 = dssp2 =='G'
zeros_and_ones2 = 1* (boolean2_1+boolean2_2)
helix_prob2 = zeros_and_ones2.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir6/parts/helix_500ns_nowat6.xtc'
#top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp3 = md.compute_dssp(traj)

boolean3_1 = dssp3 =='H'
boolean3_2 = dssp3 =='G'
zeros_and_ones3 = 1* (boolean3_1+boolean3_2)
helix_prob3 = zeros_and_ones3.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir19/parts/helix_500ns_nowat19.xtc'
#top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp4 = md.compute_dssp(traj)

boolean4_1 = dssp4 =='H'
boolean4_2 = dssp4 =='G'
zeros_and_ones4 = 1* (boolean4_1+boolean4_2)
helix_prob4 = zeros_and_ones4.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/parts/fahelix_500ns_nowat.xtc'
top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp5 = md.compute_dssp(traj)

boolean5_1 = dssp5 =='H'
boolean5_2 = dssp5 =='G'
zeros_and_ones5 = 1* (boolean5_1+boolean5_2)
helix_prob5 = zeros_and_ones5.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir2/parts/fahelix_500ns_nowat2.xtc'
top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp6 = md.compute_dssp(traj)

boolean6_1 = dssp6 =='H'
boolean6_2 = dssp6 =='G'
zeros_and_ones6 = 1* (boolean6_1+boolean6_2)
helix_prob6 = zeros_and_ones6.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir6/parts/fahelix_500ns_nowat6.xtc'
top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp7 = md.compute_dssp(traj)

boolean7_1 = dssp7 =='H'
boolean7_2 = dssp7 =='G'
zeros_and_ones7 = 1* (boolean7_1+boolean7_2)
helix_prob7 = zeros_and_ones7.sum(axis=0)

traj_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir19/parts/fahelix_500ns_nowat19.xtc'
top_file = '/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp8 = md.compute_dssp(traj)

boolean8_1 = dssp8 =='H'
boolean8_2 = dssp8 =='G'
zeros_and_ones8 = 1* (boolean8_1+boolean8_2)
helix_prob8 = zeros_and_ones8.sum(axis=0)




#plt.clf()
fig, ax = plt.subplots(2,figsize=(9,10),dpi=300)
#plt.figure(figsize=(9,10), dpi=300)
plt.subplots_adjust(hspace=.0)
ax[0].plot(helix_prob1/50000,color='blue',label='290K')
ax[0].plot(helix_prob2/50000,color='green',label='300K')
ax[0].plot(helix_prob3/50000,color='orange',label='320K')
ax[0].plot(helix_prob4/50000,color='red',label='405K')
ax[0].legend(loc="upper left",fontsize=22)
ax[1].set_xlabel('Residue Index',fontsize=24)
ax[1].set_xticks(np.arange(0, 174, 20))
ax[1].tick_params(labelsize=22)
ax[0].tick_params(labelsize=22)
#plt.yticks(np.arange(0,0.8,0.2),fontsize=20)
ax[0].set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax[0].set_xticks([])
#ax.yaxis.set_ticks([0.0,0.8,1.0])
ax[0].yaxis.set_ticks([0.0,0.2,0.4,0.6,0.8,1.0])
#ax.yaxis.set_labels([0.0,0.8,1.0])
#ax.yaxis.set_labels([0.0,0.2,0.4,0.6,0.8])
#plt.tick_params(which='both', width=2)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)

ax[0].set_xlim([0,174])
ax[0].set_ylim([0,1.0])
ax[0].set_ylabel('Helical Propensity',fontsize=24)

ax[1].set_ylim([0,1.0])
ax[1].set_xlim([0,174])
ax[1].set_yticks([0.0,0.2,0.4,0.6,0.8])
ax[1].yaxis.set_ticks([0.0,0.2,0.4,0.6,0.8])
ax[0].yaxis.set_label_coords(-.1, 0)
ax[1].plot(helix_prob5/50000,color='blue')
ax[1].plot(helix_prob6/50000,color='green')
ax[1].plot(helix_prob7/50000,color='orange')
ax[1].plot(helix_prob8/50000,color='red')

ax[0].axvspan(63,73,color='yellow',alpha=0.5, lw=0)
ax[0].axvspan(95.5,96.5,color='red',alpha=0.5,lw=0)
ax[0].axvspan(113,120,color='deepskyblue',alpha=0.5, lw=0)
ax[0].axvspan(137,142,color='deepskyblue',alpha=0.5, lw=0)
ax[0].axvspan(150,155,color='deepskyblue',alpha=0.5, lw=0)

ax[1].axvspan(63,73,color='yellow',alpha=0.5, lw=0)
ax[1].axvspan(95.5,96.5,color='red',alpha=0.5,lw=0)
ax[1].axvspan(113,120,color='deepskyblue',alpha=0.5, lw=0)
ax[1].axvspan(137,142,color='deepskyblue',alpha=0.5, lw=0)
ax[1].axvspan(150,155,color='deepskyblue',alpha=0.5, lw=0)


plt.savefig('test.png',bbox_inches="tight")

