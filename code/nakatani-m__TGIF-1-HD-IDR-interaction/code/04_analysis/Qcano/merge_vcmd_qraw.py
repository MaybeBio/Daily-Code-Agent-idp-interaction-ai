#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import copy
import kkmm_vcmd

def opt_parse():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_qcano',
                 help="q_cano file in the previous run")
    p.add_option('-o', dest='fn_out',
                 help="filename for output")
    p.add_option('--o-qraw', dest='fn_out_qraw',
                 help="filename for output")
    p.add_option('--o-qraw-cnt', dest='fn_out_qraw_cnt',
                 help="filename for output")
    p.add_option('--i-qraw', dest='fn_qraw',
                 action="append",
                 help="q_raw files")
    p.add_option('--i-qraw-list', dest='fn_qraw_list',
                 help="q_raw files")
    p.add_option('--i-weight', dest='fn_weight',
                 help="weight for each trajectory")
    p.add_option('--p-count', dest='pseudo_count',
                 type="int", default = 0,
                 help="pseudo_count")
    p.add_option('--symmetrize', dest='symmetrize',
                 action="store_true",
                 help="Symmetrize the counts")
    p.add_option('--out-log-param', dest='out_log_param',
                 action="store_true",
                 help="Output paramters in log value")
    opts, args = p.parse_args()
    print("----------------------------")
    p.print_help()
    print("----------------------------")
    return opts, args

def read_fnlist(fn):
    f = open(fn)
    fnlist = []
    for line in f:
        tmp = line.strip().split()[0]
        fnlist.append(tmp)
    f.close()
    return fnlist

def read_weight(fn):
    f = open(fn)
    weights = []
    for line in f:
        tmp = line.strip().split()[0]
        weights.append(float(tmp))
    f.close()
    return weights

def _main():
    opts, args = opt_parse()
    vc = kkmm_vcmd.VcMDConf()
    fn_list = []
    if opts.fn_qraw:
        fn_list = opts.fn_qraw
    if opts.fn_qraw_list:
        fn_list.extend(read_fnlist(opts.fn_qraw_list))
    
    trj_weight = []
    for i in range(len(fn_list)):
        trj_weight.append(1.0)
    if opts.fn_weight:
        trj_weight = read_weight(opts.fn_weight)


    vc.read_params(fn_list[0], False)
    tmp = vc.sum_params()
    vc.scale_params(trj_weight[0])
    print("%30s : samples %15.1f (%15.1f) : weight %6.4f"%(fn_list[0], tmp, vc.sum_params(), trj_weight[0]))

    for i, fn_qraw in enumerate(fn_list[1:]):
        vc_sub = kkmm_vcmd.VcMDConf()
        vc_sub.read_params(fn_qraw, False)
        tmp = vc_sub.sum_params()
        vc_sub.scale_params(trj_weight[i+1])        
        print("%30s : samples %15.1f (%15.1f) : weight %6.4f"%(fn_qraw, tmp, vc_sub.sum_params(), trj_weight[i+1]))
        vc.add_params(vc_sub)

    if opts.symmetrize:
       vc.symmetrize()

    if opts.fn_out_qraw_cnt:
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out_qraw_cnt).write(vc, param_type=0)

    if opts.fn_out_qraw:
        vc_tmp = copy.deepcopy(vc)
        vc_tmp.normalize_params()
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out_qraw).write(vc_tmp, param_type=0)

    vc.add_const(opts.pseudo_count)
    vc.normalize_params()

    vc.set_default_param()
    vc.statistics()

    vc_prev = kkmm_vcmd.VcMDConf()
    print("dbg0709c")
    vc_prev.read_params(opts.fn_qcano, False)
    if opts.symmetrize:       vc_prev.symmetrize()
    vc.multiply_params(vc_prev)
    vc.normalize_params()
    #vc.pop_zero_vs()

    vc.set_default_param()

    if opts.fn_out:
        param_mode = 0
        if opts.out_log_param:
            param_mode = 1
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out).write(vc, param_type=0, param_mode=param_mode)



if __name__ == "__main__":
    _main()
