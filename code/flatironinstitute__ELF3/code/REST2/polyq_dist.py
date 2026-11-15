####Plots the polyQ distances for 7Q by temperature, HCG ensembles
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

pq290=np.loadtxt('/mnt/home/jlindsay/ceph/elf3prd/hcg/analysis/dist/7q/290K/polyq/pqdist.dat')
pq300=np.loadtxt('/mnt/home/jlindsay/ceph/elf3prd/hcg/analysis/dist/7q/300K/polyq/pqdist.dat')
pq320=np.loadtxt('/mnt/home/jlindsay/ceph/elf3prd/hcg/analysis/dist/7q/320K/polyq/pqdist.dat')
pq415=np.loadtxt('/mnt/home/jlindsay/ceph/elf3prd/hcg/analysis/dist/7q/415K/polyq/pqdist.dat')

n_bins=50


##How to make the bars fill their bins? (not be skinny)
plt.figure(figsize=(8.5,6.25), dpi=300)

plt.hist([pq290],bins=n_bins,label=['290K'],edgecolor='black')
plt.hist([pq300],bins=n_bins,label=['300K'],edgecolor='black')
plt.hist([pq320],bins=n_bins,label=['320K'],edgecolor='black')
plt.hist([pq415],bins=n_bins,label=['415K'],edgecolor='black')

plt.legend(loc='upper left',fontsize=20)
plt.xlabel('C-alpha Distance (nm)', fontsize=20)
plt.ylabel('# Structures', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('polyQ_dist.png')
