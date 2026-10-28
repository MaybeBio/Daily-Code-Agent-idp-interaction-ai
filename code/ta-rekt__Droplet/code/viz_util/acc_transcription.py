# data_dir = "data/20260309/two_enhancer_search_nobulk_-0.3"
# data_dirs = ["data/20260310/two_enhancer_search_bulk_0.005_nuc",
#              "data/20260310/two_enhancer_search_bulk_0.010_nuc",
#              "data/20260310/two_enhancer_search_bulk_0.015_nuc"]
# data_dirs = ["data/20260403/two_enhancer_search_bulk_0.0",
#             "data/20260403/two_enhancer_search_bulk_0.010",
#             "data/20260403/two_enhancer_search_bulk_0.010_nuc"]
data_dirs = ["data/20260413/proximal_promoter_transc_d_50_c_0.003",
             "data/20260413/proximal_promoter_transc_d_50_c_0.004",
             "data/20260413/proximal_promoter_transc_d_50_c_0.005",
             "data/20260413/proximal_promoter_transc_d_50_c_0.006",
             "data/20260413/proximal_promoter_transc_d_50_c_0.007",
             "data/20260413/proximal_promoter_transc_d_50_c_0.008",
             "data/20260413/proximal_promoter_transc_d_50_c_0.009",
             "data/20260413/proximal_promoter_transc_d_50_c_0.010",
             "data/20260413/proximal_promoter_transc_d_50_c_0.011",
             "data/20260413/proximal_promoter_transc_d_50_c_0.012",
             "data/20260413/proximal_promoter_transc_d_50_c_0.013",
             "data/20260413/proximal_promoter_transc_d_50_c_0.014",
             "data/20260413/proximal_promoter_transc_d_50_c_0.015",
             "data/20260413/proximal_promoter_transc_d_50_c_0.016",
             "data/20260413/proximal_promoter_transc_d_50_c_0.017",
             "data/20260413/proximal_promoter_transc_d_50_c_0.018"]

# dists = [0]*len(data_dirs)
inds = range(1, 17)
xons, xoffs, xbs, xunbs, dists = [], [], [], [], []
dists, transcs, binds = [], [], []
n1 = 10
n = 5000

for k, data_dir in enumerate(data_dirs):
    xb, xunb = np.array([]),np.array([])
    xon, xoff = np.array([]),np.array([])
    dist, transc, bind = np.array([]), np.array([]), np.array([])
    xons.append(xon)
    xoffs.append(xoff)
    xbs.append(xb)
    xunbs.append(xunb)
    dists.append(dist)
    transcs.append(transc)
    binds.append(bind)

    for alpha_index in inds:
        path = data_dir+'/output_cpu{}.json'.format(alpha_index)

        with open(path, 'r') as f:
            res_ = json.load(f)
            res = res_[str(alpha_index)][str(1)]

        xons[k] = np.concatenate([xons[k], res["on_times"][n1:]])
        xoffs[k] = np.concatenate([xoffs[k], res["off_times"][n1:]])
        xbs[k] = np.concatenate([xbs[k], res["bound_times"][n1:]])
        xunbs[k] = np.concatenate([xunbs[k], res["unbound_times"][n1:]])
        dists[k] = np.concatenate([dists[k], res['pattern_distances'][n:]])
        transcs[k] = np.concatenate([transcs[k], res['transcription_activity'][n:]])
        binds[k] = np.concatenate([binds[k], res['binding_activity'][n:]])