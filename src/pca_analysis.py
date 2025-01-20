from numpy import array, zeros, diag, sign, argmax, zeros_like, copy, argsort, outer, triu_indices_from
from numpy.linalg import svd, eigh
import numpy as np
from .seq_io import trim_msa, encode_align, encode_align_pairs, pair_freq
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse import csr_matrix as sparsify


def compute_covariance_mat(seql, wei=None):
    bin_seq = encode_align(seq_l)
    N, Lq = bin_seq.shape
    if wei is None:
        wei = np.ones((1, N))/N
    sbin = sparsify(bin_seq)
    f_i = wei.dot(np.array(sbin.todense()))[0]
    f_ij = (sbin.T.dot(diags(wei[0], 0)).dot(sbin)).todense()
    return f_ij, f_i


def compute_pca(cov_mat):
    eigenvalues, eigenvectors = eigh(cov_mat)
    sorted_indices = argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    principal_components = eigenvectors
    explained_variance = eigenvalues
    return principal_components, explained_variance
