from svd_src.seq_io import read_fasta, encode_align, pair_freq
from svd_src.pca_analysis import compute_covariance_mat, compute_pca
from svd_src.dca_tools import compute_mf_dca, apply_apc
from svd_src.dca_tools import get_dca_scores, apply_apc, k_dca_scores
from svd_src.struct import stability_struct, read_struct, pdb_contacts
from svd_src.struct import paired_positions, get_matrix_from_db, get_matrix_pairs
from tools.mc_nested import compute_bpp, compute_ens_defect
from svd_src.data import map_rfam_pdb
import numpy as np
import matplotlib.pyplot as plt
from RNA import fold_compound
import argparse
from os.path import basename


def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('rfam')
    parser.add_argument('--struct', action="store_true")
    parser.add_argument('--vrna', action="store_true")
    parser.add_argument('-r', '--reg', type=float, default=0.1)
    return parser.parse_args()


def main():
    args = parse_arguments()
    fam = basename(args.rfam)[:-4]
    msa = read_fasta(f"./data/coconet/RNA_DATASET/MSA/{fam}.faclean")

    seq_l = [seq for seq in msa.values()]
    nb_seq = len(seq_l)
    L = len(seq_l[0])
    Lq = L * 5
    bin_seq = encode_align(seq_l).astype(np.int16)

    if args.struct:
        output_f = f"out/mf_dca_{fam}_3.dca"
        if args.vrna:
            struct = read_struct(f"./dca_input/vrna_{fam}.str")
        else:
            struct = read_struct(f"./data/coconet/RNA_DATASET/secstruct/secstruct_{map_rfam_pdb[fam]}.txt")
        nrj_l = np.array([stability_struct(seq, struct) for seq in seq_l])
        nrj_l = (nrj_l - nrj_l.mean())/nrj_l.std()
        temp = 1
        nrj_w = np.exp((1./temp) * nrj_l)
        if nrj_w.sum() > 1e-12:
            nrj_w = nrj_w/nrj_w.sum()
        nb_seq = np.exp(-(nrj_w * np.log(nrj_w)).sum())
    else:
        output_f = f"out/mf_dca_{fam}_0.dca"
        nrj_w = None

    dca_rewei, pf_r = compute_mf_dca(bin_seq, wei_nrj=nrj_w, reg=0.01)
    dca_rewei_apc = apply_apc(dca_rewei)

    with open(output_f, "w") as out:
        out.write(f"# EFF SEQ: {nb_seq}\n")
        for i in range(L):
            for j in range(i+1, L):
                out.write("{} {} {}\n".format(i, j, dca_rewei[i, j]))

if __name__ == '__main__':
    main()
