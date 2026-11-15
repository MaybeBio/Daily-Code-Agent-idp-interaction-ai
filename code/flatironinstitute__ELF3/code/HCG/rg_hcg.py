####This script creates a plot of all temperature conditions of 7Q HCG

import mdtraj as md
import matplotlib.pyplot as plt

##Load in ensembles as trajectories

t290=md.load('/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/290K/7q_290K_all.xtc',top='/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/290K/1/elf3prd_7q_290K_set1_conf0.pdb')
t300=md.load('/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/300K/7q_300K_all.xtc',top='/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/300K/1/elf3prd_300K_set1_conf0.pdb')
t320=md.load('/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/320K/7q_320K_all.xtc',top='/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/320K/1/elf3prd_7q_320K_set1_conf0.pdb')
t415=md.load('/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/415K_pseudo/415Kp_all.xtc',top='/mnt/home/jlindsay/ceph/elf3prd/hcg/pdbs/7q/415K_pseudo/1/elfprd_415Kp_set1_conf0.pdb')

t290rg=md.compute_rg(t290)
t300rg=md.compute_rg(t300)
t320rg=md.compute_rg(t320)
t415rg=md.compute_rg(t415)

n_bins=50

#fig, axs = plt.subplots(1,2,sharex=True, sharey=True, tight_layout=True)
#axs[0].hist(rgwt,bins=n_bins)
#axs[1].hist(rgmut,bins=n_bins)
#plt.savefig('fig.png')

plt.figure(figsize=(8.5,6), dpi=300)

plt.hist([t290rg],n_bins,label=['290K'],edgecolor='black')
plt.hist([t300rg],n_bins,label=['300K'],edgecolor='black')
plt.hist([t320rg],n_bins,label=['320K'],edgecolor='black')
plt.hist([t415rg],n_bins,label=['415K'],edgecolor='black')

plt.legend(loc='upper right',fontsize=20)
plt.xlabel('Gyradius (nm)', fontsize=20)
plt.ylabel('# Structures', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('fig.png')

