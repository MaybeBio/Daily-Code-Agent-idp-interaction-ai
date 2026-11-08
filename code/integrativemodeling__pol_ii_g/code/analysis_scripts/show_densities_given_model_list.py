import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import os,sys
import glob

#precision_cutoff = float(sys.argv[1]) # precision threshold for calculating resolution of the MRC
path="Clusone.0/"

densityc={"GDOWN11":['GDOWN1'],"Nt":[(1,94,'GDOWN1')],"Ct":[(300,335,'GDOWN1')],"C2":[(216,300,'GDOWN1')],"N2":[(1,65,'GDOWN1')]} 


gmd1 = IMP.pmi.analysis.GetModelDensity(custom_ranges=densityc, resolution=12.2, voxel=3.0)

for rmfile in glob.glob("%s/*.rmf3" % path):
    model = IMP.Model()
    inf = RMF.open_rmf_file_read_only(rmfile)
    h = IMP.rmf.create_hierarchies(inf, model)[0]
    #if int(rmfile.split("/")[-1].split("_")[0]) < 250:
    gmd1.add_subunits_density(h)
    #else:
        #gmd2.add_subunits_density(h)

gmd1.write_mrc(path)
