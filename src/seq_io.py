"""Read files
"""

from numpy import triu_indices_from, diag_indices_from, array, argmax, int16

NUC = ["A", "C", "G", "U", "-"]
NUC_P = [(a, b) for a in ['A', 'C', 'G', 'U', '-']
         for b in ['A', 'C', 'G', 'U', '-']]
NUC_B = [
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1],
]

def num_to_seq(num_seq):
    "convert num to seq"
    return "".join([NUC[el] for el in num_seq])


def seq_to_index(seq):
    res = []
    for el in seq:
        res += [NUC.index(el) if el in NUC else NUC.index("-")]
    return res


def msa_to_index(msa):
    res_msa = []
    for seq in msa:
        res_msa += [seq_to_index(seq)]
    return array(res_msa)


def encode_seq(seq, neg=False):
    res = []
    for el in seq:
        if neg:
            res += NNUC_B[NUC.index(el)] if el in NUC else NNUC_B[NUC.index("-")]
        else:
            res += NUC_B[NUC.index(el)] if el in NUC else NUC_B[NUC.index("-")]
    return res


def decode_seq(bmsa, neg=False):
    msa = []
    bsize, seq_len_5 = bmsa.shape
    if neg:
        bmsa = bmsa.reshape(-1, seq_len_5//4, 4)
    else:
        bmsa = bmsa.reshape(-1, seq_len_5//5, 5)
    for batch in bmsa:
        msa += ["".join([NUC[n] for n in argmax(batch, axis=-1)])]
    return msa


def encode_align(msa, neg=False):
    res_msa = []
    for seq in msa:
        res_msa += [encode_seq(seq, neg)]
    return array(res_msa)


def encode_conf_full(seq, st_pairs):
    "One-hot encoding of paired positions"
    conf = [0 for _ in range(len(st_pairs)*25)]
    for ii, (pi, pj) in enumerate(st_pairs):
        if (seq[pi], seq[pj]) in NUC_P:
            conf[ii*25 + NUC_P.index((seq[pi], seq[pj]))] = 1
    return conf


def encode_align_pairs(msa, struct):
    "struct = [(pi, pj)]"
    res_msa = []
    for seq in msa:
        res_msa += [encode_conf_full(seq, struct)]
    return array(res_msa)



def read_fasta(infile):
    results = {}
    for l in open(infile):
        if l.startswith(">"):
            name = l.strip()[1:]
            results[name] = ""
        else:
            results[name] += l.strip()
    return results


def read_seq(infile):
    results = {}
    for si, l in enumerate(open(infile)):
        results[si] = l.strip()
    return results


def trim_msa(msa, ref_seq):
    pos = [i for i, si in enumerate(ref_seq) if si != "-"]
    if type(msa) is dict:
        new_msa = {}
        for n, seq in msa.items():
            new_msa[n] = "".join([seq[i] for i in pos])
    else:
        new_msa = []
        for seq in msa:
            new_msa += ["".join([seq[i] for i in pos])]
    return new_msa, "".join([ref_seq[i] for i in pos])


def pair_freq(bmsa, no_gap:bool=False, diag_zero=False):
    if no_gap:
        nb_seq, seq_sp = bmsa.shape
        bmsa = bmsa.reshape(-1, seq_sp//5, 5)[:, :, :-1].reshape(-1, 4*(seq_sp//5))
    Lq = bmsa.shape[1]
    nb_seq, seq_sp = bmsa.shape
    res = (bmsa[:, :, None] * bmsa[:, None, :]).mean(axis=0).reshape(Lq, Lq)
    return res


def single_freq(bmsa, no_gap:bool=False):
    if no_gap:
        nb_seq, seq_sp = bmsa.shape
        bmsa = bmsa.reshape(-1, seq_sp//5, 5)[:, :, :-1].reshape(-1, 4*(seq_sp//5))
    return bmsa.mean(axis=0)


def reweight_msa(bmsa, cor_thres=0.8):
    cor = bmsa @ bmsa.T/(bmsa.shape[1]//5)
    iwei = array([sum(cor[p, :] > cor_thres) for p in range(bmsa.shape[0])])
    return bmsa/iwei.reshape(bmsa.shape[0], 1)
