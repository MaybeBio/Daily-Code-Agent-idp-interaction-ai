import mdtraj as md
import matplotlib.pyplot as plt

top='/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/helix_nowat.pdb'
traj='/mnt/home/jlindsay/ceph/elf3prd/rest/helix/prod/dir0/ongoing/helix_nowat200ns.xtc'
wt=md.load(traj,top=top)
rgwt=md.compute_rg(wt)
N_points=len(rgwt)
top='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/fahelix_nowat.pdb'
traj='/mnt/home/jlindsay/ceph/elf3prd/rest/fa_helix/prod/dir0/ongoing/fahelix_nowat200ns.xtc'
mut=md.load(traj,top=top)
rgmut=md.compute_rg(mut)
n_bins=100

#fig, axs = plt.subplots(1,2,sharex=True, sharey=True, tight_layout=True)
#axs[0].hist(rgwt,bins=n_bins)
#axs[1].hist(rgmut,bins=n_bins)
#plt.savefig('fig.png')

plt.figure(figsize=(8.5,6), dpi=300)

plt.hist([rgwt],n_bins,label=['WT'],edgecolor='black',alpha=0.5)
plt.hist([rgmut],n_bins,label=['F96A'], edgecolor='black',alpha=0.5)

plt.legend(loc='upper right',fontsize=20)
plt.xlabel('Gyradius (nm)', fontsize=20)
plt.ylabel('# Structures', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('fig.png')
