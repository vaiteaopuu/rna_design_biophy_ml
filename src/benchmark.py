from svd_src.dca_tools import get_dca_scores, apply_apc, k_dca_scores, read_alpha
from svd_src.pca_analysis import compute_covariance_mat, compute_pca
from svd_src.seq_io import read_fasta, encode_align, pair_freq
from random import sample
from RNA import fold_compound
import numpy as np
from scipy.stats import gaussian_kde

def read_plm_cocnet(rfam):
    infile = f"./data/coconet/RAW_COEV_DATA_ALL/RAW_DATA_DATASET/CoCoNet_RAW_DATA_FiveFold_CrossVal/TwoFilterMatrices/trial_1/fold_1/{rfam}__2x7x7.txt"
    res = {}
    pos = set()
    for l in open(infile):
        if not l.startswith("#") and len(l.strip()) > 0:
            pi, pj, score = l.strip().split()
            res[(int(pi), int(pj))] = float(score)
            pos.add(int(pi))
            pos.add(int(pj))
    mat = np.zeros((max(pos), max(pos)))
    for (pi, pj), sc in res.items():
        mat[pi-1, pj-1] = sc
        mat[pj-1, pi-1] = sc
    return mat

def read_plm_cocnet(rfam):
    infile = f"./data/coconet/RAW_COEV_DATA_ALL/RAW_DATA_DATASET/plmDCA_RAW_DATA_DATASET/{rfam}_plmDCA.txt"
    res = {}
    pos = set()
    for l in open(infile):
        if not l.startswith("#"):
            pi, pj, score = l.strip().split()
            res[(int(pi), int(pj))] = float(score)
            pos.add(int(pi))
            pos.add(int(pj))
    mat = np.zeros((max(pos), max(pos)))
    for (pi, pj), sc in res.items():
        mat[pi-1, pj-1] = sc
        mat[pj-1, pi-1] = sc
    return mat


def read_bl_cocnet(rfam):
    infile = f"./data/coconet/RAW_COEV_DATA_ALL/RAW_DATA_DATASET/BL_RAW_DATA_DATASET/{rfam}_BL.txt"
    res = {}
    pos = set()
    for l in open(infile):
        if not l.startswith("#"):
            pi, pj, score = l.strip().split()
            res[(int(pi), int(pj))] = float(score)
            pos.add(int(pi))
            pos.add(int(pj))
    mat = np.zeros((max(pos), max(pos)))
    for (pi, pj), sc in res.items():
        mat[pi-1, pj-1] = sc
        mat[pj-1, pi-1] = sc
    return mat


def read_stat(infile):
    results = []
    for l in open(infile):
        if not l.startswith("#"):
            pi, pj, ai, aj, prob, ccov = l.strip().split()
            if pi != pj:
                results += [float(prob)]
    return results


def read_data(dca_file):
    dca_mat = apply_apc(get_dca_scores(dca_file))
    alpha = read_alpha(dca_file)
    return dca_mat, alpha


def one_test(dca_mat, struct_mat=None, mask=None):
    if struct_mat is not None:
        dca_mat[struct_mat==1] = dca_mat.max()

    if mask is not None:
        dca_mat = dca_mat * (1. - mask)
    return dca_mat

def get_wt_data(fam):
    wt = []
    wt_d = {}
    for i in range(2):
        tmp_seq = [l.strip() for l in open(f"./cluster_in_2/input/{fam}/clust_{i}.seq")]
        wt += tmp_seq
        wt_d[i] = encode_align(tmp_seq)
    len_seq = len(wt[0])
    bin_seq_r = encode_align(wt)
    bin_seq = bin_seq_r.astype(float) - bin_seq_r.mean(axis=0)
    pca, eigv = compute_pca(bin_seq.T @ bin_seq)
    return wt_d, bin_seq_r, pca, eigv

def get_data_R0(ffam, bin_mean, pca):
    des_d = {}
    struct_d = {}
    des_db = {}
    for clust in [0, 1]:
        name = f"{clust}_R"
        path = f"./cluster_in_2/{ffam}/clust_{clust}_0.seq_Tclust_{clust}"
        des = [l.strip() for l in open(path)]
        des = sample(des, 500)
        struct_d[name] = [fold_compound(s).mfe()[0] for s in des]
        bin_dseq = encode_align(des).astype(float) - bin_mean
        des_d[name] = encode_align(des).astype(float)
        des_db[name] = bin_dseq @ pca
    return des_d, des_db

def compute_pvalue(gen_proj, wt_proj, num_dim):
    kde = gaussian_kde(wt_proj[:, :num_dim].T)
    random_samples = kde.resample(10000).T  # sample from WT
    likelihoods = kde(random_samples.T)     # compute likelihood of the samples
    new_likelihoods = kde(gen_proj[:, :num_dim].T)  # compute the likelihood of the generated seq
    cdf_values = np.mean(likelihoods < new_likelihoods[:, None], axis=1)
    p_values_two_tailed = 2 * np.minimum(cdf_values, 1 - cdf_values)
    return p_values_two_tailed
