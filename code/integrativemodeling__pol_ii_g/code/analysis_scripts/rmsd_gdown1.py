import IMP
import RMF
import IMP.atom
import IMP.rmf
import os,sys,math,numpy
import pickle
import IMP.pmi
import IMP.pmi.analysis
import glob
import IMP.pmi
import IMP.pmi.analysis
from scipy.spatial.distance import cdist
import numpy as np


tgt="Gdown1"
pklFile = tgt+'.cluster_0.pkl' 

f1=open('Gdown1.txt', 'w+')

def get_dic_from_two_lists(keys, values):
    return { keys[i] : values[i] for i in range(len(keys)) }

conform=[]
num=0
conform_min=[]
rmsdc=["GDOWN1"]


for file in glob.glob("GSM23/RMF/*.rmf3"):
   print >>f1, file, num
   print num,file
   partcoord= []
   labels=[]
   m = IMP.Model()
   inf = RMF.open_rmf_file_read_only(file)
   h = IMP.rmf.create_hierarchies(inf, m)[0]
   particle2s = IMP.core.get_leaves(h)
   IMP.rmf.load_frame(inf, 0)
   for state in h.get_children():
        for component in state.get_children():
                if component.get_name() in rmsdc:
                        for leaf in IMP.core.get_leaves(component):
                                p1=IMP.core.XYZ(leaf.get_particle())
                                partcoord.append(p1.get_coordinates())
  # except IOError:
#	pass
   conform.append(partcoord)
   num += 1
   partcoord = []


distmat=numpy.zeros((num+1,num+1))
avg_dist=[]
for i in range(num):
    for j in range(i+1,num):
         dist=IMP.algebra.get_rmsd(conform[i],conform[j])
         print dist,i,j
         distmat[i][j]=dist
         avg_dist.append(dist)

print sum(avg_dist),len(avg_dist),sum(avg_dist)/len(avg_dist)
pickle.dump(distmat,open(pklFile,"wb"))                                                                                                                                                                         
