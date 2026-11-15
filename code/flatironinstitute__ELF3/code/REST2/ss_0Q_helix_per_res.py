import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

traj_file='/mnt/home/jlindsay/ceph/elf3prd/rest/0Q/dir0/elf3_0Q_nowat0.xtc'
top_file='/mnt/home/jlindsay/ceph/elf3prd/rest/0Q/dir0/traj_nowat.pdb'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp1 = md.compute_dssp(traj)

boolean1_1 = dssp1 =='H'
boolean1_2 = dssp1 =='G'
zeros_and_ones1 = 1* (boolean1_1+boolean1_2)
helix_prob1 = zeros_and_ones1.sum(axis=0)

traj_file='/mnt/home/jlindsay/ceph/elf3prd/rest/0Q/dir2/elf3_0Q_nowat2.xtc'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp2 = md.compute_dssp(traj)

boolean2_1 = dssp2 =='H'
boolean2_2 = dssp2 =='G'
zeros_and_ones2 = 1* (boolean2_1+boolean2_2)
helix_prob2 = zeros_and_ones2.sum(axis=0)

traj_file='/mnt/home/jlindsay/ceph/elf3prd/rest/0Q/dir4/elf3_0Q_nowat4.xtc'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp3 = md.compute_dssp(traj)

boolean3_1 = dssp3 =='H'
boolean3_2 = dssp3 =='G'
zeros_and_ones3 = 1* (boolean3_1+boolean3_2)
helix_prob3 = zeros_and_ones3.sum(axis=0)

traj_file='/mnt/home/jlindsay/ceph/elf3prd/rest/0Q/dir11/elf3_0Q_nowat11.xtc'
traj = md.load(traj_file,top=top_file)
topology = traj.topology
dssp4 = md.compute_dssp(traj)

boolean4_1 = dssp4 =='H'
boolean4_2 = dssp4 =='G'
zeros_and_ones4 = 1* (boolean4_1+boolean4_2)
helix_prob4 = zeros_and_ones4.sum(axis=0)

#plt.clf()
fig, ax = plt.subplots()
plt.figure(figsize=(9,5), dpi=300)
plt.plot(helix_prob1/50001,color='black',label='290K')
plt.plot(helix_prob2/50001,color='red',label='300K')
plt.plot(helix_prob3/50001,color='green',label='320K')
plt.plot(helix_prob4/50001,color='blue',label='405K')
plt.legend(loc="upper left",fontsize=22)
plt.xlabel('Residue Index',fontsize=22)
plt.xticks(np.arange(0, 174, 20),fontsize=20)
#plt.yticks(np.arange(0,0.8,0.2),fontsize=20)
plt.yticks([0.0,0.2,0.4,0.6,0.8 ],fontsize=20)
#ax.yaxis.set_ticks([0.0,0.8,1.0])
#ax.yaxis.set_ticks([0.0,0.2,0.4,0.6,0.8])
#ax.yaxis.set_labels([0.0,0.8,1.0])
#ax.yaxis.set_labels([0.0,0.2,0.4,0.6,0.8])
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.axvspan(63,73,color='yellow',alpha=0.5, lw=0)
#plt.axvspan(95.5,96.5,color='red',alpha=0.5,lw=0)
plt.axvspan(130,134,color='deepskyblue',alpha=0.5, lw=0)
plt.axvspan(143,147,color='deepskyblue',alpha=0.5, lw=0)
plt.xlim([0,174])
plt.ylim([0,0.8])
plt.axvline(x=113, color = 'deepskyblue', linestyle='--')
plt.ylabel('Helical Propensity',fontsize=22)
plt.savefig('temp.png',bbox_inches="tight")

