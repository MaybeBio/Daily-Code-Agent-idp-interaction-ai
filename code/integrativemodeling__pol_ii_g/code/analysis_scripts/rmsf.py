import time
import resource
import numpy as np
import sys, os, glob

import scipy as sp
from scipy import spatial

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import os,sys
 
from math import sqrt
import IMP.rmf
import RMF

conform = []
labels=[]
weight = []
num = 0


tmxl=open("RMSF.csv",'w')


for file in glob.glob("Clusone.0/*.rmf3"):
   print num
   partcoord = []
   labels=[]
   m = IMP.Model()
   inf = RMF.open_rmf_file_read_only(file)
   h = IMP.rmf.create_hierarchies(inf, m)[0]
   particle2s = IMP.core.get_leaves(h)
   IMP.rmf.load_frame(inf, 0)
   for state in h.get_children():
      for component in state.get_children():
        for leaf in IMP.core.get_leaves(component):
            p=IMP.core.XYZ(leaf.get_particle())
            labels.append([p.get_particle_index(),component.get_name()])
            partcoord.append(p.get_coordinates())   
   conform.append(partcoord)
   partcoord = []
   num = num + 1



std =np.array(conform).std(0)*np.array(conform).std(0)                                                            
part=np.size(std)/3
for i in range(0,part):
    Xs = np.sum(std[i,0])                                
    Ys = np.sum(std[i,1])                                
    Zs = np.sum(std[i,2])
    print >> tmxl, labels[i][0],labels[i][1],sqrt((Xs+Ys+Zs))

