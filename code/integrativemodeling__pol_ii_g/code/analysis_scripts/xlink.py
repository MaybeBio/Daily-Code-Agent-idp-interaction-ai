#!/usr/bin/env python

import os,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import string
import numpy


class Violations(object):

    def __init__(self, xlinkfile,threshold):

        self.violation_threshold  = threshold 
        self.violation_counts = {}  # dictionary with a key per restraint

        xf = open(xlinkfile)
        for ln in xf.readlines(): 
            if ln.startswith('PROTEIN1'):
                continue
            (res1,pos1,res2,pos2) =ln.strip().split(',')[0:4]
        
            self.violation_counts[(res1,pos1,res2,pos2)]=0 
        xf.close()

    def get_xlink_distances(self, lndict,xlinkType):
        # get the minimum distances for each xlink'ed pair in a single model

        minimum_distances = {}
        dist_over_threshold = {}
        mindist_copies = {}

        for k in self.violation_counts:
            minimum_distances[k] = 1000.0 # some big number
            dist_over_threshold[k] = 0.0 
            mindist_copies[k] = ""
        
        for rst in lndict:
            print (rst)
            items = rst.split('|')
            if items[0]!='CrossLinkingMassSpectrometryRestraint_Distance_' or items[1]!=xlinkType:
                continue
            (res1,pos1,res2,pos2) = items[3:7]
                     
            if '.' in res1:
                resonly1 = res1.split('.')[0]  # separate the residue name from copy number
            else:
                resonly1 = res1
            if '.' in res2:
                resonly2 = res2.split('.')[0]   # separate the residue name from copy number
            else:
                resonly2=res2
        
            modeldistance = float(lndict[rst])
            #print res1,pos1,res2,pos2,resonly1,resonly2

            if modeldistance<minimum_distances[(resonly1,pos1,resonly2,pos2)]:
                minimum_distances[(resonly1,pos1,resonly2,pos2)]=modeldistance
                mindist_copies[(resonly1,pos1,resonly2,pos2)]=res1+"|"+pos1+"|"+res2+"|"+pos2 #write in the appropriate way
                #print mindist_copies[(resonly1,pos1,resonly2,pos2)],minimum_distances[(resonly1,pos1,resonly2,pos2)]                    
            
        for k in minimum_distances:
            if minimum_distances[k]<=self.violation_threshold:
                dist_over_threshold[k] = 0.0                
        
            else:
                dist_over_threshold[k] = minimum_distances[k]-self.violation_threshold
	   	    
        return minimum_distances,mindist_copies,dist_over_threshold

statFile = sys.argv[1] # stat file for the top model
xlinkFile =sys.argv[2]
xlinkType = sys.argv[3]
xlinkCutoff = float(sys.argv[4])

Analysis = Violations(xlinkFile,xlinkCutoff) 

sf=open(statFile)
for i,ln in enumerate(sf.readlines()):
        distances,mindist_copies,dist_over_threshold = Analysis.get_xlink_distances(eval(ln.strip('\n' )),xlinkType)
        print (i)
        print (ln)



# 2. get the distances for each xlink, output to a file
#sf = open(statFile)
##########ln has all the headings
#ln = sf.readlines()[0]
#print (eval(ln.strip('\n' ))
#distances,mindist_copies,dist_over_threshold = Analysis.get_xlink_distances(eval(ln.strip('\n' )),xlinkType)
#sf.close()

#tmxl = open('top_model_distances_'+xlinkType+'_cutoff-'+str(xlinkCutoff)+'.txt','w')

#print >>tmxl,"#Xlink_with_copy","Prot1","Res1","Prot2","Res2","Xlink_distance","Xlink_violation"
#for xlink in distances:
#    print >>tmxl, xlinkType+"|"+mindist_copies[xlink],xlink[0],xlink[1],xlink[2],xlink[3],"%.3f" %(distances[xlink]),"%.3f" %(dist_over_threshold[xlink])
