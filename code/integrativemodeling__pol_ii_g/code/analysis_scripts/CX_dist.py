import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.em
import IMP.core
import IMP.atom
import IMP.rmf
import RMF
import sys
import os

import glob
import numpy as np
from numpy.linalg import svd

###Reading CX file and saving it
###################################################################
###################################################################
sel1=[];sel2=[];pos1=[];pos2=[]
filename=sys.argv[1]
sf=open(filename,'r')
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split(',')
    sel1.append(line[0])
    sel2.append(line[2])
    pos1.append(line[1])
    pos2.append(line[3])
size=len(sel1)
print size
###############################################################


#####Write information to file
#################################################################
#########################
tmxl = open('Gdowndistall_clus0.txt','w')

#######################################
indir='/Clusone.0/RMF'

for root, dirs,filenames in os.walk(indir):
    for rmf_name in filenames:
        print (rmf_name)
        m = IMP.Model()
        f = RMF.open_rmf_file_read_only(os.path.join(root,rmf_name))
        bps = IMP.rmf.create_hierarchies(f, m)[0]
        ps = IMP.core.get_leaves(bps)
        IMP.rmf.load_frame(f, 0)
        distances=[];violations=[]
        for i in range(size):
            p1 = IMP.atom.Selection(bps,
                                    molecule =sel1[i],
                                    residue_index=int(pos1[i])).get_selected_particles()
            
	    c1 = IMP.core.XYZR(p1[0]);X1=float(c1.get_x());Y1=float(c1.get_y());Z1=float(c1.get_z())
            p2 = IMP.atom.Selection(bps,
                                    molecule =sel2[i],
                                    residue_index=int(pos2[i])).get_selected_particles()
	   
            c2 = IMP.core.XYZ(p2[0]);X2=float(c2.get_x());Y2=float(c2.get_y());Z2=float(c2.get_z())
            dist=((X1-X2)**2 + (Y1-Y2)**2 +(Z1-Z2)**2)**0.5
	    print dist
	    print >> tmxl, dist

