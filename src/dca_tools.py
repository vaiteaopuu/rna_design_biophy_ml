"""This module contains function to analyze DCA results
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_matrix as sparsify


def read_alpha(infile):
    for l in open(infile):
        if l.startswith("# ALPHA"):
            return float(l.strip().split()[-1])


def read_dca_bl_scores(infile):
    "read the output of bl-dca"
    results = {}
    for line in open(infile):
        nrj, pi, pj = line.split()
        results[(int(pi), int(pj))] = float(nrj)
    return results


def read_dca_score(infile):
    "read the output of ALF-RNA"
    results = {}
    for line in open(infile):
        if not line.startswith("#"):
            posi_, posj_, nuc_i, nuc_j, nrj_ = line.strip().split()
            posi, posj = int(posi_), int(posj_)
            nrj = float(nrj_)
            if "-" not in [nuc_i, nuc_j]:
                if (posi, posj) not in results:
                    results[(posi, posj)] = {(nuc_i, nuc_j): nrj}
                else:
                    results[(posi, posj)][(nuc_i, nuc_j)] = nrj
    return results

def comp_dca_scores(seq, dca):
    "compute dca score"
    tot = 0.0
    for i, ai in enumerate(seq):
        for j, aj in enumerate(seq[i:], start=i):
            tot += dca[(i, j)][(ai, aj)]
    return tot

def compute_frobenius_norm(couplings, positions):
    "from the couplings compute the frobenius norm"
    frob = {}
    frob_i = {}
    results = {}
    tot = 0.0

    for i, j in couplings:
        frob[(i, j)] = sum(nrj**2 for (nuc_a, nuc_b), nrj in couplings[(i, j)].items())**0.5
        tot += sum(nrj**2 for (nuc_a, nuc_b), nrj in couplings[(i, j)].items())**.5

    for pos in positions:
        frob_i[pos] = sum(el**2 for pair, el in frob.items() if pos in pair)**0.5

    for i, j in frob:
        results[(i, j)] = frob[(i, j)] - frob_i[i]*frob_i[j]/tot
    return results


def compute_frobenius_norm_mat(couplings, positions):
    "from the couplings compute the frobenius norm"
    len_seq = len(positions)
    frob = np.zeros((len_seq, len_seq))

    for i, j in couplings:
        if abs(i - j) > 3:
            frob[i, j] = sum(nrj**2 for (nuc_a, nuc_b), nrj in couplings[(i, j)].items())**0.5

    return frob


def get_dca_scores(parm_files):
    dca = read_dca_score(parm_files)
    positions = list(set([pi for pi, pj in dca]))
    positions.sort()
    return compute_frobenius_norm_mat(dca, positions)


def compute_cov(bin_seq, wei=None):
    N, Lq = bin_seq.shape
    if wei is None:
        wei = np.ones((1, N))/N
    sbin = sparsify(bin_seq)
    f_i = wei.dot(np.array(sbin.todense()))[0]
    f_ij = (sbin.T.dot(diags(wei[0], 0)).dot(sbin)).todense()
    return f_ij, f_i


def compute_mf_dca(msa, q=5, reg=0.01, wei_nrj=None, bsize=100):
    N, Lq = msa.shape  # N: number of sequences, Lq: L * q
    L = Lq // q  # Number of positions
    if wei_nrj is not None:
        if len(wei_nrj.shape) != 2 or wei_nrj.shape[0] != 1:
            wei_nrj = wei_nrj[None, :]

    f_ij, f_i = compute_cov(msa, wei_nrj)

    # Centered covariance matrix
    C = f_ij - np.outer(f_i, f_i)

    # Regularization
    C += reg * np.eye(Lq)

    # Invert the covariance matrix
    J = -np.linalg.inv(C)

    # Calculate DCA scores
    dca_scores = np.zeros((L, L))
    for i in range(L):
        for j in range(i+1, L):
            # Calculate Frobenius norm for the (i, j) block in J
            submatrix = J[i*q:(i+1)*q, j*q:(j+1)*q]
            dca_scores[i, j] = np.linalg.norm(submatrix, 'fro')
            dca_scores[j, i] = dca_scores[i, j]  # Symmetric matrix

    return dca_scores, f_ij - np.outer(f_i, f_i)


def apply_apc(C):
    # Row and column means
    row_mean = np.mean(C, axis=1)
    col_mean = np.mean(C, axis=0)

    # Overall mean
    overall_mean = np.mean(C)

    # APC correction
    C_apc = C - np.outer(row_mean, col_mean) / overall_mean

    return C_apc


def k_dca_scores(dca_score, k):
    norm = np.argpartition(dca_score.flatten(), -k)[-k:]
    i, j = np.unravel_index(norm, dca_score.shape)
    return [el for el in zip(i, j)]


def seq_reweight(bin_seq, max_seqid=0.8):
    seq_len = bin_seq.shape[1]//21
    sbin = sparsify(bin_seq)
    sim_mat = sbin.dot(sbin.T) / seq_len
    seqw = np.array(1 / (sim_mat > max_seqid).sum(axis=0))
    return seqw
