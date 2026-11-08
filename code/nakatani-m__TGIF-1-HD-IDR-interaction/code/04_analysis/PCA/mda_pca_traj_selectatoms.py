#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import sys
import os
import re
import copy
from optparse import OptionParser

import kkstruct
import kkpdb
import kkmdconf
from MDAnalysis import Universe
from MDAnalysis.analysis import distances
import numpy as np

def define_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_mdconf',
                 help="file name for MDConf file")
    p.add_option('--trr', dest='fn_trrs',
                 action="append",
                 help="file names for trajectories")
    p.add_option('--ref', dest='fn_ref',
                 help="file names for reference pdb")
    p.add_option('-o', dest='fn_out',
                 help="file name for output")

    p.add_option('--o-eig', dest='fn_out_eig',
                 help="file name for output")
    p.add_option('--o-crd', dest='fn_out_crd',
                 help="file name for output coordinates")

    p.add_option('--select', dest='selection',
                 default = "name *",
                 help="atom selection string")
    p.add_option('--dim', dest='n_dim',
                 default=-1, type="int",
                 help="dimension for output")
    
    p.add_option('--i-eig', dest='fn_in_eig',
                 help="input eigen vectors")

    p.add_option('--begin', dest='frame_begin',
                 type = "int", default=0,
                 help="the first frame number")
    p.add_option('--end', dest='frame_end',
                 type = "int", default=-1,
                help="the last frame number")
    p.add_option('--range', dest='frame_range',
                 action="append",
                help="the frame range: i-j")

    p.add_option('--i-prob', dest='fn_prob',
                help="reweighted probability file")
    p.add_option('--col-prob', dest='col_prob',
                 type="int",
                 help="The column-ID for weights in --i-prob file")
    p.add_option('--non-log-prob', dest='flg_non_log_prob',
                 action="store_true",
                 help="The input probabilities is non-log values.")
    p.add_option('--min-prob', dest='min_prob',
                 type = "float", default=0.0,
                help="The minimum threshold of probability")

    p.add_option('--o-prob', dest='fn_o_prob',
                help="Output prob file")

    p.add_option('--coordinate-type', dest='coordinate_type',
                 type="choice",
                 choices=["cartesian","distance", "bitvector"],
                 default="cartesian",
                 help="")
    p.add_option('--atom-pairs', dest='atom_pairs',
                 action="append",
                 help='"i:j" means pair i and j, "i-j:k-l" means all-vs-all pairs between i to j and k to l')
    p.add_option('--fn-atom-pairs', dest='fn_atom_pairs',
                 help="")
    p.add_option('--interval-snap', dest='interval_snap',
                 type="int", default=1,
                 help="")

    p.add_option('--skip-pca', dest='flg_skip_pca',
                 action="store_true",
                 help="Skip PCA")

    p.add_option('--vec-col', dest='vec_col',
                 type="int",
                 help="The column id for bit vector")
    p.add_option('--vec-header', dest='vec_header',
                 action="store_true",
                 help="THe header line exists in the bit-vector file")

    opts, args = p.parse_args()
    print("----------------------------")
    p.print_help()
    print("----------------------------")
    return opts, args


def _main():
    opts, args = define_options() 
    fn_init = opts.fn_ref
    list_trr = opts.fn_trrs
    if opts.fn_mdconf:
        mdconf = kkmdconf.MDConfReader(opts.fn_mdconf).read_conf()
        fn_init = mdconf.gro_init_structure
        list_trr = mdconf.gro_traj_trr

    weights = []
    if opts.fn_prob:
        print("read_weight: min:" + str(opts.min_prob))
        weights = read_weight(opts.fn_prob, opts.fn_o_prob, 
                              opts.min_prob,
                              not opts.flg_non_log_prob,
                              opts.col_prob)
        print(len(weights))

    range_info = []
    if opts.frame_range:
        range_info = parse_frame_range(opts.frame_range)
    else:
        range_info.append((opts.frame_begin, opts.frame_end))

    print("read_trajectory")
    crd_vector = []
    if opts.coordinate_type == "cartesian":
        crd_vector = read_trajectory(fn_init, list_trr, opts.selection, 
                                     weights, 
                                     opts.min_prob,
                                     range_info,
                                     opts.interval_snap)
    elif opts.coordinate_type == "distance":
        atom_pairs_pre = []
        if opts.fn_atom_pairs:
            atom_pairs_pre.extend(read_atom_pairs(opts.fn_atom_pairs))
        if len(opts.atom_pairs) > 0:
            atom_pairs_pre.extend(opts.atom_pairs)
        print("atom_pairs")
        print(atom_pairs_pre)
        atom_pairs = parse_atom_pair(atom_pairs_pre)
        print(atom_pairs)
        crd_vector = read_trajectory_dist(fn_init, list_trr, atom_pairs,
                                          weights, 
                                          opts.min_prob,
                                          range_info,
                                          opts.interval_snap)
    elif opts.coordinate_type == "bitvector":
        crd_vector = read_bit_vector(list_trr, weights, opts.min_prob, range_info, 
                                     opts.vec_col, opts.vec_header,
                                     opts.interval_snap)
        
        
    print("crd_vector")
    print(crd_vector)
    if opts.fn_out_crd:
        write_crd(opts.fn_out_crd, crd_vector)
        
    mu = None
    if opts.fn_in_eig:
        eigval, eigvec, mu = read_eig(opts.fn_in_eig)
    elif not opts.flg_skip_pca:
        idx, eigval, eigvec, mu = calc_pca(crd_vector)
        if opts.fn_out_eig:
            output_eig(opts.fn_out_eig, idx, eigval, eigvec, mu,  opts.selection)

    print("EIGEN VALUES")
    print(eigval)
    print("EIGEN VALUES (RATIO)")
    print(eigval/sum(eigval))
    print("EIGEN VECTORS")
    print(eigvec)
    output_pca_traj(opts.fn_out, eigval, eigvec, mu, crd_vector, n_dim=opts.n_dim)

def read_weight(fn_prob, fn_o_prob, min_prob, flg_exp=True, col_weight=5):
    """
    Read the weighting coefficients of each snapshot
    'fn_prob' is the file name of space-separated table,
      which records the weights as 'col_weight'-th column.

    ## Obsoleted 
    #Returns the new crd_vector, which is a subset of the argument crd_vector,
    #  by extracting snapshots with the weight >= min_prob.
    #The subset of fn_prob with >= min_prob is written on for fn_o_prob.
    """
    
    ## Read the weights 
    weight_buf = []
    f_prob = open(fn_prob)
    for line in f_prob:
        if line[0] == "#": continue
        terms = line.strip().split()
        weight_buf.append(float(terms[col_weight]))
    f_prob.close()
    weight_ar = np.array(weight_buf)
    if flg_exp: weight_ar = np.exp(weight_ar)
    max_weight = weight_ar.max()
    #assert(crd_vector.shape[0] == len(weight_ar))

    #new_crd_vector = []
    #for i, crd in enumerate(crd_vector):
    #    if weight_ar[i]/max_weight >= min_prob:
    #        tmp = []
    #        for val in crd: tmp.append(val)
    #        new_crd_vector.append(tmp)
    #ret_crd_vector = np.array(new_crd_vector)
    #print ret_crd_vector[0]

    ## output prob subset
    fo_prob = None
    if fn_o_prob: fo_prob = open(fn_o_prob,"w")
    f_prob = open(fn_prob)
    n_prob = 0
    for line in f_prob:
        if line[0] == "#":
            if fo_prob: fo_prob.write(line)
            continue
        terms = line.strip().split()
        p = float(terms[col_weight])
        if flg_exp: p = np.exp(p)
        if p / max_weight >= min_prob:
            ## print terms[col_weight] + " exp:" + str(np.exp(float(terms[col_weight]))) + " " + str(np.exp(float(terms[col_weight]))/max_weight) + " " + str(min_prob)
            if fo_prob: fo_prob.write(line)
            n_prob += 1
    f_prob.close()
    if fo_prob: fo_prob.close()
    print("frames p>=" + str(min_prob) + " : " + str(n_prob))

    #print "# Snapshots : " + str(ret_crd_vector.shape[0]) + " / " + str(crd_vector.shape[0])
    
    #return ret_crd_vector
    return weight_ar/max_weight

def read_trajectory(fn_init, trr_list, selection,
                    weights, min_prob,
                    range_info, frame_step=1):
    crd_vector = []
    univ = None
    if len(trr_list) == 0:
        univ = Universe(fn_init)
        for i_frame, ts in enumerate(univ.trajectory):
            if i_frame%frame_step != 0: continue
            if len(weights) > 0 and weights[i_frame] < min_prob:
                continue
            atoms = univ.select_atoms(selection)
            tmp_vec = []
            for crd in atoms.positions:
                tmp_vec.extend(crd)
            crd_vector.append(tmp_vec)
    else:
        i_frame = -1
        for fn_trr in trr_list:
            univ = Universe(fn_init, fn_trr)

            for i, ts in enumerate(univ.trajectory):
                i_frame += 1
                flg_read = False
                for rg in range_info:
                    if rg[0] <= i_frame and \
                       (rg[1] < 0 or i_frame <= rg[1]): flg_read = True
                if not flg_read: continue
                if i%frame_step != 0: continue
                if len(weights) > 0 and weights[i_frame] < min_prob:
                    continue
                atoms = univ.select_atoms(selection)
                tmp_vec = []
                for crd in atoms.positions:
                    tmp_vec.extend(crd)
                crd_vector.append(tmp_vec)
    #print(crd_vector.shape)
    return np.array(crd_vector)

def read_trajectory_dist(fn_init, trr_list, atom_pairs,
                         weights, min_prob,
                         range_info, frame_step=1):
    """
    Read the trajectory and calculating interatomic distances.
    The argument 'groups' is list of pairs of atoms: [(i,j), (k,l), ...]
    """
    crd_vector = None
    univ = None
    if len(trr_list) == 0:
        tmp_vector = []
        univ = Universe(fn_init)
        for i, ts in enumerate(univ.trajectory):
            if i%frame_step != 0: continue
            if len(weights) > 0 and weights[i_frame] < min_prob:
                continue
            tmp_vec = []
            #for crd in univ.atoms.coordinates():
            atoms = univ.atoms.positions
            for atpair in atom_pairs:
                crd0 = atoms[atpair[0]]
                crd1 = atoms[atpair[1]]
                dist = np.sqrt(((crd0-crd1)**2).sum())
                tmp_vec.append(dist)
            print(tmp_vec)
            tmp_vector.append(tmp_vec)
        crd_vector = np.array(tmp_vector)
    else:
        # counting the number of frames
        n_frames = 0
            
        i_frame = -1
        for fn_trr in trr_list:
            univ = Universe(fn_init, fn_trr)
            for i, ts in enumerate(univ.trajectory):
                i_frame += 1
                if len(weights) > 0 and weights[i_frame] < min_prob: continue
                if i%frame_step != 0: continue
                flg_read = True
                for rg in range_info:
                    if rg[0] > i_frame or \
                       (rg[1] >= 0 and i_frame >= rg[1]): flg_read = False
                if flg_read: n_frames += 1
        print("n_frames : " + str(n_frames))

        crd_vector = np.zeros((n_frames, len(atom_pairs)), dtype=np.float32)
        i_frame = -1
        i_frame_dat = -1
        for fn_trr in trr_list:
            univ = Universe(fn_init, fn_trr)
            
            for i, ts in enumerate(univ.trajectory):
                i_frame += 1
                flg_read = True
                for rg in range_info:
                    if rg[0] > i_frame or \
                       (rg[1] >= 0 and i_frame >= rg[1]): flg_read = False
                if not flg_read: continue
                #if frame_begin > i_frame: continue
                #if frame_end > 0 and frame_end < i_frame: continue
                if i%frame_step != 0: continue
                if len(weights) > 0 and weights[i_frame] < min_prob:
                    continue
                #atoms = univ.selectAtoms(selection)
                i_frame_dat += 1
                atoms = univ.atoms.positions
                #print "dbg " + str(i_frame_dat)
                #print atoms
                tmp_vec = []
                for i, atpair in enumerate(atom_pairs):
                    dist = 0
                    try:
                        dist = np.sqrt(((atoms[atpair[0]]-atoms[atpair[1]])**2).sum())
                    except:
                        print("EXCEPT:")
                        print(atpair)
                        print(atoms)
                        sys.exit(1)
                        
                    crd_vector[i_frame_dat, i] = dist
                #print crd_vector[i_frame_dat]
                    #tmp_vec.append(dist)
                #print tmp_vec
                #crd_vector.append(tmp_vec)
                #print "vec " + str(len(tmp_vec)) + " / " + str(len(crd_vector))
    #print "length crd_vector: " + str(len(crd_vector))
    #return np.array(crd_vector)
    return crd_vector

def read_bit_vector(list_vec, weights, min_prob, range_info, interval_snap, vec_col, vec_header, frame_step=1):
    n_frames = 0
    i_frame = -1
    for fn_vec in list_vec:
        f = open(fn_vec)
        if vec_header:
            header = f.readline()
        for line in f:
            i_frame += 1
            if len(weights) > 0 and weights[i_frame] < min_prob: continue
            if weights[i_frame] == 0 : continue
            if i%frame_step != 0: continue
            flg_read = True
            for rg in range_info:
                if rg[0] > i_frame or \
                   (rg[1] >= 0 and i_frame >= rg[1]): flg_read = False
            if not flg_read: continue
            vc = np.array(list(line.strip().split()[vec_col]), dtype=np.int8)
            tmp_vector.append(tmp_vector)
        f.close()
        crd_vector = np.array(tmp_vector)
    return crd_vector

def read_atom_pairs(fn_atom_pairs):
    """
    Read input file specifying atom pairs for distance calculations
    """
    ap = []
    f = open(fn_atom_pairs)
    for line in f:
        ap.append(line.strip())
    f.close()
    return ap

def parse_atom_pair(atom_pairs_pre):
    """
    Parsing atom pair definition strings.
    "0:1" means a atom pair 0 and 1. This is converted into a tuple (0,1).
    "0-1:2-3" means all-vs-all pairs of (0,1) and (2,3).
    This is converted into a list of tuples: [(0,2),(0,3),(1,2),(1,3)]
    """
    atom_pairs = []
    #print atom_pairs_pre
    for pair in atom_pairs_pre:
        re_v = re.compile("^([\d\-\,]+)\:([\d\-\,]+)$")
        re_m = re_v.match(pair)
        if not re_m:
            sys.stderr.write("Invalid atom pair definition (:): \n")
            sys.stderr.write(pair + "\n")
            continue
        grp_pre = [re_m.group(1), re_m.group(2)]
        #print grp_pre
        grp = [[],[]]
        for i, grp_pre_str in enumerate(grp_pre):
            terms = grp_pre_str.split(",")
            for j, term in enumerate(terms):
                re_v2 = re.compile("^(\d+)(\-(\d+)){0,1}$")
                re_m2 = re_v2.match(term)
                if not re_m2:
                    sys.stderr.write("Invalid atom definition (-):\n")
                    sys.stderr.write(grp_pre_str+"\n")
                    sys.stderr.write(terms+"\n")
                    break
                if re_m2.group(3):
                    grp[i].extend(range(int(re_m2.group(1)), int(re_m2.group(3))+1))
                elif re_m2.group(1):
                    grp[i].append(int(re_m2.group(1)))
        #print grp[0]
        #print grp[1]
        for a1 in grp[0]:
            #print "a1 " + str(a1)
            for a2 in grp[1]:
                #print "a2 " + str(a2)
                a1n = a1
                if a1n > a2:
                    tmp = a2; a2 = a1n; a1n = tmp
                elif a1n==a2: continue
                #print (a1n,a2)
                if not (a1n, a2) in atom_pairs:
                    atom_pairs.append((a1n,a2))
    return atom_pairs

def parse_frame_range(frame_ranges):
    """
    Parsing frame-range definition strings.
    Arg:
      frame_ranges : a list of strings
        An element specifies the range of frames to be read
        e.g., "10-20"
    Return:
      range_info : a list of tuples
        An element is two-integer tuple, defining range of frames
    """

    range_info = []
    for rng in frame_ranges:
        re_v = re.compile("^(\d+)\-(\d+)$")
        re_m = re_v.match(rng)
        if not re_m:
            sys.stderr.write("Invalid atom pair definition:")
            sys.stderr.write(pair)
            continue
        range_info.append((int(re_m.group(1)), int(re_m.group(2))))
    return range_info

def calc_pca(crd_vector):
    crd_m = np.mean(crd_vector, axis=0)
    #crd_sd = np.sqrt(np.var(crd_vector, axis=0))
    #crd_z = (crd_vector-crd_m)/crd_sd
    crd_n = (crd_vector-crd_m).T
    cov = np.dot(crd_n, crd_n.T) / float(crd_vector.shape[0])
    #l, v = np.linalg.eig( np.dot(crd_z.T, crd_z) )
    l, v = np.linalg.eig( cov )
    idx = l.argsort()
    l = l[idx][::-1]
    v = v[:,idx][:,::-1]
    return idx, l, v, crd_m

def output_eig(fn_out, idx, eigval, eigvec, mu, selection):
    #print idx
    #univ = Universe(mdconf.gro_init_structure)
    #    atoms = univ.selectAtoms(selection)
    #resids = []
    #for atom in atoms:
    #    for dim in ["x","y","z"]:
    #        resids.append(str(atom.resid) + "." + dim)
    #resids = resids[idx][::-1]
    f = open(fn_out,"w")
    line = "\t".join([str(x) for x in mu])
    f.write(line+"\n")
    line = "\t".join([str(x) for x in eigval])
    f.write(line+"\n")
    
    #f.write("\t".join(resids)+"\n")
    for eig in eigvec:
        line = "\t".join([str(x) for x in eig])
        f.write(line+"\n")        
    f.close()
    return 0

def read_eig(fn_in):
    f = open(fn_in)
    mu = [float(x) for x in re.compile("\s+").split(f.readline().strip())]
    eigval = [float(x) for x in re.compile("\s+").split(f.readline().strip())]
    eigvec = []
    for line in f:
        terms = [float(x) for x in re.compile("\s+").split(line.strip())]
        eigvec.append(terms)
    f.close()
    return np.array(eigval), np.array(eigvec), np.array(mu)

def output_pca_traj(fn_out, eigval, eigvec, mu, crd_vector, n_dim = -1):
    #crd_m = np.mean(crd_vector, axis=0)
    #crd_sd = np.sqrt(np.var(crd_vector, axis=0))
    #crd_z = (crd_vector-crd_m)/crd_sd
    
    f = open(fn_out,"w")
    mtx = np.dot(eigvec.T, (crd_vector - mu).T).T
    print("mtx.shape: " + str(mtx.shape[0]) + " - " + str(mtx.shape[1]))
    if n_dim < 0:
        n_dim = len(mtx[0])
    for i, row in enumerate(mtx):
        line = str(i) + "\t"
        line += "\t".join([str(x) for x in row[0:n_dim]])
        f.write(line+"\n")
    #for i,ts in enumerate(crd_z):
    #    line = str(i) + "\t"
    #    line += "\t".join([str(x) for x in np.dot(ts, eigvec)[0:n_dim]])
    #    f.write(line+"\n")
    f.close()
    return 0

def write_crd(fn_out_crd, crd_vector):
    f = open(fn_out_crd, "w")
    for i_frm in range(crd_vector.shape[0]):
        line = "%10d"%(i_frm+1)
        #f.write(line+"\n")
        #line = "%16.9f"%(1.0)
        for i_dim in range(crd_vector.shape[1]):
            line += "%16.9f"%crd_vector[i_frm, i_dim]
        f.write(line+"\n")
    f.close()
    return

if __name__=="__main__":
    _main()
