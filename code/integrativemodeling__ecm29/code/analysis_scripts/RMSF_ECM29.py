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
weight = []
num = 0

path=str(sys.argv[1])
for file in glob.glob("%s/*.rmf3" % path):
   print num
   partcoord = []
   m = IMP.Model()
   inf = RMF.open_rmf_file_read_only(file)
   h = IMP.rmf.create_hierarchies(inf, m)[0]
   particle2s = IMP.core.get_leaves(h)
   IMP.rmf.load_frame(inf, 0)
   for state in h.get_children():
      for component in state.get_children():
         if "ecm29" in component.get_name():
            for leaf in IMP.core.get_leaves(component):
               if len(IMP.atom.Fragment(leaf).get_residue_indexes()) == 0:
                  p=IMP.core.XYZ(leaf.get_particle())
                  partcoord.append(p.get_coordinates())
                  if num == 0:
                     weight.append(1.0)
               else:
                  p=IMP.core.XYZ(leaf.get_particle())
                  partcoord.append(p.get_coordinates())
                  if num == 0:
                     weight.append(float(len(IMP.atom.Fragment(leaf).get_residue_indexes())))
   
   conform.append(partcoord)
   partcoord = []
   num = num + 1

std =np.array(conform).std(0)*np.array(conform).std(0)                                                            
Xs = np.sum(np.dot(std[:,0], np.array(weight))) / np.sum(np.array(weight))                                
Ys = np.sum(np.dot(std[:,1], np.array(weight))) / np.sum(np.array(weight))                                
Zs = np.sum(np.dot(std[:,2], np.array(weight))) / np.sum(np.array(weight))

print sqrt((Xs+Ys+Zs))

print ""
print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
print ""
