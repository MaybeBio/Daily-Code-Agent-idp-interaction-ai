import os,sys,random,numpy,math
import pickle
import scipy as sp
from scipy import spatial
import scipy.stats

def get_run_identity(idfile):
    # whether a run belongs to run1 or run2
    run1_models=[]
    run2_models=[]

    with open(idfile, 'r') as f:
        for line in f:
            mod = line.split("/")[-1]
            if int(mod.split("_")[0]) < 250:
                run1_models.append(int(mod.strip("\n").split(" ")[1]))
            else:
                run2_models.append(int(mod.strip("\n").split(" ")[1]))

    return run1_models, run2_models

def get_random_models(all_models,num_models):
	final=[]

	while len(final)<num_models:
		sel=random.choice(all_models) 

                if sel not in final:
			final.append(sel)
                        all_models.remove(sel) 
	return final

def get_cutoffs_list(distmat,gridSize,numModels):

	maxdist=0.0
	mindist=1000000.0
	for i in range(numModels-1):
    		for j in range(i+1,numModels):
        		if distmat[i][j]>maxdist:
            			maxdist=distmat[i][j]

        		if distmat[i][j]<mindist:
           			 mindist=distmat[i][j]

	cutoffs=numpy.arange(mindist,maxdist,gridSize) # or maxdist/2.0, 5.0 s

	return cutoffs

def get_subset_distance_matrix(distmat_full, models_subset):

	distmat_subset=numpy.zeros((len(models_subset),len(models_subset)))

	for i in range(len(models_subset)-1):
   		for j in range(i+1,len(models_subset)):
			model_index_i=models_subset[i]
			model_index_j=models_subset[j]
                        
            		if model_index_i<model_index_j:
				distmat_subset[i][j]=distmat_full[model_index_i][model_index_j]
			else:
				distmat_subset[i][j]=distmat_full[model_index_j][model_index_i]

	return distmat_subset

def precision_cluster(distmat,numModels,rmsd_cutoff):
    #STEP 2. Populate the neighbors ofa given model
    neighbors=[]
    for count in range(numModels):
        neighbors.append([count])  # model is a neighbor of itself
 
    for i in range(numModels-1):
        for j in range(i+1,numModels):
        
            if distmat[i][j]<=rmsd_cutoff: # accepted to be a neighbor
                #print i,j,distmat[i][j]

                neighbors[i].append(j)
                neighbors[j].append(i)
     
    #STEP 3. Get the weightiest cluster, and iterate
    unclustered=[]
    boolUnclustered=[]
    for i in range(numModels):
        unclustered.append(i)
        boolUnclustered.append(True)

    cluster_members=[] # list of lists : one list per cluster
    cluster_centers=[]

    while len(unclustered)>0:
        # get cluster with maximum weight
        max_neighbors=0
        currcenter=-1
        for eachu in unclustered:  # if multiple clusters have same maxweight this tie is broken arbitrarily! 
            if len(neighbors[eachu])>max_neighbors:
                max_neighbors=len(neighbors[eachu])
                currcenter=eachu   
   
        #form a new cluster with u and its neighbors
        cluster_centers.append(currcenter)
        cluster_members.append([n for n in neighbors[currcenter]]) 

        #update neighbors 
        for n in neighbors[currcenter]:
            #removes the neighbor from the pool
            unclustered.remove(n) #first occurence of n is removed. 
            boolUnclustered[n]=False # clustered

        for n in neighbors[currcenter]:
            for unn in neighbors[n]: #unclustered neighbor
                if not boolUnclustered[unn]:
                    continue
                neighbors[unn].remove(n)
   
    return cluster_centers,cluster_members

def get_contingency_table(num_clusters,cluster_members,models_subset,run1_models,run2_models):

	full_ctable=numpy.zeros((num_clusters,2))
		
	for ic,cluster in enumerate(cluster_members):
		for member in cluster:
			model_index=models_subset[member]

			if model_index in run1_models:
                                #print "run1",model_index
                                full_ctable[ic][0]+=1.0
			elif model_index in run2_models:
				#print "run2",model_index
                                full_ctable[ic][1]+=1.0

	## now normalize by number of models in each run
	numModelsRun1 = float(numpy.sum(full_ctable,axis=0)[0])
	numModelsRun2 = float(numpy.sum(full_ctable,axis=0)[1])

        ##print numModelsRun1, numModelsRun2
	#for i in range(num_clusters):
		
		#full_ctable[i][0]=full_ctable[i][0]*100.0/numModelsRun1
			
		#full_ctable[i][1]=full_ctable[i][1]*100.0/numModelsRun2

	reduced_ctable=[]
	
	retained_clusters=[]

        num_models_run1_10percent=0.10*numModelsRun1

        num_models_run2_10percent=0.10*numModelsRun2

	for i in range(num_clusters):
		if full_ctable[i][0]<=num_models_run1_10percent and full_ctable[i][1]<=num_models_run2_10percent:
			continue
		reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
		retained_clusters.append(i)

	return numpy.array(reduced_ctable),retained_clusters

def test_sampling_convergence(contingency_table,total_num_models):

    if len(contingency_table)==0:
        return 0.0,1.0

    ct = numpy.transpose(contingency_table) 
    
    [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)
    
    if dof==0.0:
            cramersv=0.0 #converged, one single cluster
    else:
            cramersv=math.sqrt(chisquare/float(total_num_models))
   
             
    return(pvalue,cramersv)

def percent_ensemble_explained(ctable,total_num_models):
       
        if len(ctable)==0:
            return 0.0
        percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)
        
        return percent_clustered

# each run of this script gives, for a randoms subset of num_models_needed, the p-value 
tgt=sys.argv[1]
idfile=sys.argv[2] # gives the run identity for a cluster
pklFile=sys.argv[3]
total_num_models=int(sys.argv[4])

# parameters to set
gridSize=2.5

#GetModel Lists
run1_all_models,run2_all_models=get_run_identity(idfile)

#Load RMSD matrix
import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
mHandler = MatrixHandler()
mHandler.loadMatrix(pklFile)
rmsd_matrix = mHandler.getMatrix()
distmat = rmsd_matrix.get_data()  
distmat_full = sp.spatial.distance.squareform(distmat)

#GetCutOffs
cutoffs_list=get_cutoffs_list(distmat_full,gridSize,total_num_models)
print cutoffs_list

all_models=run1_all_models+run2_all_models
num_models_per_run=total_num_models/2.0
'''
pvals=[]
cvs=[]
percents=[]

for c in cutoffs_list:

    cluster_centers,cluster_members=precision_cluster(distmat_full,total_num_models,c)
    ctable,retained_clusters=get_contingency_table(len(cluster_centers),cluster_members,all_models,run1_all_models,run2_all_models)
    (pval,cramersv)=test_sampling_convergence(ctable,total_num_models)
    percent_explained= percent_ensemble_explained(ctable,total_num_models)
    
    pvals.append(pval)
    cvs.append(cramersv)
    percents.append(percent_explained)
    print c
# Now apply the rule for selecting the right precision and the pvalue/cramersv
precision_convergence=100000.0

for i in range(len(cutoffs_list)):
    print cutoffs_list[i], percents[i], pvals[i], cvs[i]
    if percents[i]>80.0 and (pvals[i]>0.05 or cvs[i]<0.10):
        if precision_convergence>cutoffs_list[i]:
            precision_convergence=cutoffs_list[i]
    else: # takes care of local minima? We want to choose the first time at which the above condition was satisfied
        precision_convergence=100000.0
'''
# Redo the clustering at the required precision of sampling convergence
#TODO could have saved time by storing this beforehand
precision_convergence = 60.00
cluster_centers,cluster_members=precision_cluster(distmat_full,total_num_models,precision_convergence)

# We need to know which clusters to ignore and which to focus on!
ctable,retained_clusters=get_contingency_table(len(cluster_centers),cluster_members,all_models,run1_all_models,run2_all_models)
(pval,cramersv)=test_sampling_convergence(ctable,total_num_models)
percent_explained= percent_ensemble_explained(ctable,total_num_models)
print precision_convergence, pval, cramersv, percent_explained

# Output models to files
for i in range(len(retained_clusters)):
    clus=retained_clusters[i] 
    
    both_file=open('cluster.'+str(i)+'.all.txt','w')
    run1_file=open('cluster.'+str(i)+'.run1.txt','w')
    run2_file=open('cluster.'+str(i)+'.run2.txt','w')

    for mem in cluster_members[clus]: #TODO Not efficient, could have done this while doing the contingency table
        print >>both_file,mem
        if mem in run1_all_models:
            print >>run1_file, mem
        else:
            print >>run2_file, mem

    both_file.close()
    run1_file.close()
    run2_file.close()

