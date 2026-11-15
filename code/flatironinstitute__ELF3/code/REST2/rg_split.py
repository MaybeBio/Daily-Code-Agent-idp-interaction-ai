import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('talk')

top='/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/helix_nowat.pdb'
traj='/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/parts/helix_500ns_nowat.xtc'

wt=md.load(traj, atom_indices=[0,1477], top=top)
rgwt=md.compute_rg(wt)
N_points=len(rgwt)
mut=md.load(traj, atom_indices=[1478,2605], top=top)
rgmut=md.compute_rg(mut)
n_bins=100

top='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/parts/fahelix_500ns_nowat.xtc'
wt=md.load(traj, atom_indices=[0,1467], top=top)
rgwt2=md.compute_rg(wt)
N_points=len(rgwt2)
#top='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
#traj='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/parts/fahelix_500ns_nowat.xtc'
mut=md.load(traj, atom_indices=[1468,2595], top=top)
rgmut2=md.compute_rg(mut)
n_bins=100

fig, axs = plt.subplots(2,1,sharex=True, tight_layout=True)
plt.subplots_adjust(wspace=0)
axs[0].hist(rgwt,bins=n_bins,label="Res 1-100",color="black")
axs[0].hist(rgmut,bins=n_bins,label="Res 101-173",color="gray", alpha=0.5)
axs[1].hist(rgwt2,bins=n_bins, label="Res 1-100",color="black")
axs[1].hist(rgmut2,bins=n_bins, label="Res 101-173",color="gray",alpha=0.8)
#plt.savefig('fig.png')

#plt.figure(figsize=(8.5,6), dpi=300)
#plt.hist([rgwt],n_bins,label=['Resid 1-100'],edgecolor='black',alpha=0.5)
#plt.hist([rgmut],n_bins,label=['Resid 101-173'], edgecolor='black',alpha=0.5)
axs[0].legend(loc='upper right')
axs[1].set_xlabel('Radius of gyration (nm)')
axs[0].set_ylabel('# Structures', fontsize=20)
axs[0].text(x=5.5,y=700,s="WT 290K")
axs[1].text(x=5.5,y=450,s="F96A 290K")
axs[0].set_yticks([0,1000,2000,3000])
axs[1].set_yticks([0,1000,2000,3000])
axs[0].set_ylim([0,4000])
axs[1].set_ylim([0,4000])
#axs[1].set_xticks(fontsize=16)
#axs[0].set_yticks(fontsize=16)
plt.savefig('test.png',bbox_inches='tight')
