from Bio.PDB import PDBParser
from RNA import fold_compound, bp_distance
import numpy as np


def paired_positions(sequence):
    "give the paired positions"
    # save open bracket in piles
    pile_reg, pile_pk = [], []
    pairs = []

    try:
        for i, sstruc in enumerate(sequence):
            if sstruc in "(":
                pile_reg += [i]
            elif sstruc == "[":
                pile_pk += [i]
            elif sstruc == ")":
                pairs += [(pile_reg.pop(), i)]
            elif sstruc == "]":
                pairs += [(pile_pk.pop(), i)]
    except:
        ""
    return pairs


def get_matrix_pairs(pair_list, len_seq):
    mat = np.zeros((len_seq, len_seq))
    for i, j in pair_list:
        mat[i, j] = 1.
    return mat


def get_matrix_from_db(db_struct):
    len_seq = len(db_struct)
    mat = np.zeros((len_seq, len_seq))
    for i, j in paired_positions(db_struct):
        mat[i, j] = 1.
    return mat


def stability_struct(seq, struct):
    "compute secondary structure stability"
    return fold_compound(seq).eval_structure(struct)


def read_struct(infile):
    "read from coconet dataset"
    struct = ""
    for l in open(infile):
        l.replace(" ", "")
        if l.startswith("SECSTRUCT"):
            struct += l.strip().split(":")[1].strip()
    return struct


# XXX 3D structure
def get_heavy_atoms(residue):
    """
    Extracts heavy atoms (non-hydrogen) from a given residue.
    """
    return [atom for atom in residue if atom.element != "H"]


def remove_sec_3D(pair_sec, pair_tri):
    tmp = []
    for pi, pj in pair_tri:
        if (pi, pj) not in pair_sec and\
           (pi-1, pj) not in pair_sec and\
           (pi-2, pj) not in pair_sec and\
           (pi, pj+1) not in pair_sec and\
           (pi, pj+2) not in pair_sec:
            tmp += [(pi, pj)]
    return tmp


def compute_residue_contacts(residues, threshold=5.0):
    """
    Compute contacts between residues based on heavy atoms.
    """
    contacts = []
    num_residues = len(residues)

    for i in range(num_residues):
        for j in range(i + 3, num_residues):
            res1_atoms = get_heavy_atoms(residues[i])
            res2_atoms = get_heavy_atoms(residues[j])

            # Extract coordinates of heavy atoms in numpy arrays
            coords1 = np.array([atom.coord for atom in res1_atoms])
            coords2 = np.array([atom.coord for atom in res2_atoms])

            # Calculate pairwise distances and check for any distance below the threshold
            distances = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=-1)
            if np.any(distances <= threshold):
                contacts.append((i, j))

    return contacts


def pdb_contacts(pdb_file, threshold=5.0):
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA_structure", pdb_file)

    # Extract RNA residues
    rna_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Only include standard RNA nucleotides
                if residue.id[0] == " " and residue.resname in ["A", "C", "G", "U"]:
                    rna_residues.append(residue)

    # Compute contacts between residues
    return compute_residue_contacts(rna_residues, threshold)


def compute_bpp(seq):
    seq_comp = fold_compound("".join(seq))
    mfe_st, mfe_scr = seq_comp.pf()
    pairs_mat = np.array(seq_comp.bpp())[1:, 1:]
    return pairs_mat


def compute_ens_defect(seq, target):
    seq_comp = fold_compound("".join(seq))
    mfe_st, mfe_scr = seq_comp.pf()
    pairs_mat = np.array(seq_comp.bpp())[1:, 1:]
    return np.log(seq_comp.ensemble_defect(target))


def compute_corr(seq, target_mat):
    seq_comp = fold_compound("".join(seq))
    # seq_comp.hc_add_from_db(HC)
    mfe_st, mfe_scr = seq_comp.pf()
    pairs_mat = np.array(seq_comp.bpp())[1:, 1:]
    score = np.sum((target_mat - pairs_mat)**2)**0.5
    return score, mfe_st


def compute_corr_mod(seq, target_mat):
    seq_comp = fold_compound("".join(seq))
    # seq_comp.hc_add_from_db(HC)
    mfe_st, mfe_scr = seq_comp.pf()
    pairs_mat = np.array(seq_comp.bpp())[1:, 1:]
    score = np.sum(target_mat * pairs_mat)
    return score, mfe_st


def get_subopt(seq, delta_nrj=500):
    subopt_data = { 'counter' : 1, 'sequence' : seq , "structs": []}
    def store_subopt_result(structure, energy, data):
        if not structure == None:
            data['counter'] = data['counter'] + 1
            data['structs'] += [(structure, energy)]
    a = fold_compound(seq)
    a.subopt_cb(delta_nrj, store_subopt_result, subopt_data);
    return subopt_data["structs"]


def pairwise_dist_struct(struct_list):
    len_seq = len(struct_list[0])
    nb_struct = len(struct_list)
    dist = np.zeros((nb_struct, nb_struct))
    print(dist.shape)
    for i, si in enumerate(struct_list):
        for j, sj in enumerate(struct_list[i+1:], start=i+1):
            dist[i, j] = bp_distance(si, sj)
            dist[j, i] = dist[i, j]
    return 1-(dist - dist.min())/(dist.max() - dist.min())
