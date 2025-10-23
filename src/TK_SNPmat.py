from Bio import SeqIO
import pandas as pd

import argparse
import os

from itertools import combinations
from multiprocessing import Process, Manager

import numpy as np

from OrthoSLC import __version__, mission_spliter

p = argparse.ArgumentParser(
    prog="python3 TK_SNPmat.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)

p.add_argument("-i", "--input_path",
               metavar='',
                required=True,
                help="<dir> path/to/input/directory")
p.add_argument("-o", "--output_csv",
               metavar='',
                required=True,
                help="<dir> path/to/output/file.csv")
p.add_argument("-T", "--ID_TSV",
               metavar='',
                required=True,
                help="<txt> path/to/Step1.txt")
p.add_argument("-u", "--resource_total",
               metavar='',
                default="1",
                help="<int> total threads available")

args, extras = p.parse_known_args()

kaligned_path = args.input_path

CONV_DICT = {
    "A": [1.0, 0.0, 0.0, 0.0],
    "T": [0.0, 1.0, 0.0, 0.0],
    "C": [0.0, 0.0, 1.0, 0.0],
    "G": [0.0, 0.0, 0.0, 1.0],
    "Y": [0.0, 0.5, 0.5, 0.0],
    "R": [0.5, 0.0, 0.0, 0.5],
    "S": [0.0, 0.0, 0.5, 0.5],
    "W": [0.5, 0.5, 0.0, 0.0],
    "K": [0.0, 0.5, 0.0, 0.5],
    "M": [0.5, 0.0, 0.5, 0.0],
    "B": [0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    "D": [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0],
    "H": [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
    "V": [1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0],
    "N": [0.25, 0.25, 0.25, 0.25],
    "-": [0.0, 0.0, 0.0, 0.0],
}


def encode_sequence(sequence):
    """Convert a DNA sequence to a flattened numerical representation."""

    try:
        encoded = np.array([CONV_DICT[base.upper()] for base in sequence], dtype=np.float32)
    except KeyError as exc:
        raise ValueError(f"Unsupported base '{exc.args[0]}' encountered in sequence.") from exc

    return encoded.reshape(-1)


def pairwise_l2_distance(vectors):
    """Compute pairwise L2 distance matrix for the provided vectors."""

    if not vectors:
        return np.empty((0, 0), dtype=np.float32)

    mat = np.vstack(vectors).astype(np.float32)
    sq_norms = np.sum(mat ** 2, axis=1, keepdims=True)
    dist_sq = sq_norms - 2 * mat @ mat.T + sq_norms.T
    np.maximum(dist_sq, 0, out=dist_sq)
    return np.sqrt(dist_sq, out=dist_sq)

def distance_in_one_cluster(in_path_ls, 
                            share_ls, 
                            loop_ls, 
                            total_len):
    species_order = sorted({sid for sid_pair in loop_ls for sid in sid_pair})

    for in_path in in_path_ls:
        in_fasta = SeqIO.to_dict(SeqIO.parse(kaligned_path + '/' + in_path,
                                             'fasta'))

        spe_dict = {x[0:x.index('-')]: x for x in in_fasta.keys()}

        encoded_vectors = []
        cluster_seq_len = None
        for sid in species_order:
            record_id = spe_dict.get(sid)
            if record_id is None:
                raise KeyError(f"Sequence for species '{sid}' missing in alignment '{in_path}'.")

            seq = str(in_fasta[record_id].seq)
            if cluster_seq_len is None:
                cluster_seq_len = len(seq)

            encoded_vectors.append(encode_sequence(seq))

        distance_matrix = pairwise_l2_distance(encoded_vectors)
        species_index = {sid: idx for idx, sid in enumerate(species_order)}

        result_dict = {}

        for pair in loop_ls:
            idx_a = species_index[pair[0]]
            idx_b = species_index[pair[1]]
            result_dict[pair] = distance_matrix[idx_a, idx_b]

        total_len.append(cluster_seq_len)
        share_ls.append(result_dict)

if __name__ == '__main__':
    procc_num = int(args.resource_total)
    
    mission_ls = mission_spliter(os.listdir(kaligned_path), procc_num)
    
    manager = Manager()
    return_dict_ls = manager.list()
    len_ls = manager.list()
    
    a_fasta = SeqIO.parse(kaligned_path + '/' + os.listdir(kaligned_path)[0], 'fasta')
    a_ls = [x.id[0: x.id.index('-')] for x in a_fasta]
    a_ls.sort()
    
    par_ls = list(combinations(a_ls, 2))
    
    for submissions in mission_ls:
    
        jobs = []

        p = Process(target = distance_in_one_cluster,
                    args = (submissions, 
                            return_dict_ls,
                            par_ls,
                            len_ls
                           )

                   )
        p.start()
        jobs.append(p)


    for z in jobs:
        z.join()
    return_dict_ls = list(return_dict_ls)

#add up all the snp number
sum_snp = {}
for key in return_dict_ls[0]: #prepare a dict to save snp number
    sum_snp[key] = 0
    
for i in return_dict_ls:
    for key in sum_snp:
        sum_snp[key] = sum_snp[key] + i[key]

SLC_path_1 = args.ID_TSV

df_SLC_1 = pd.read_csv(SLC_path_1, sep = "\t"
                       , header=None
                       , index_col = 0)
dict_kalign_id = df_SLC_1[1].to_dict()

sum_snp_4_2 = {(dict_kalign_id[int(k[0])], dict_kalign_id[int(k[1])]): v for k, v in sum_snp.items()}
sum_snp_4_2 = dict(sorted(sum_snp_4_2.items()))

strain_naams = list(df_SLC_1[1])
df_mat = pd.DataFrame(0
                      , index=strain_naams
                      , columns=strain_naams
                     )
df_mat.index = strain_naams
df_mat.columns = strain_naams

cbn = list(combinations(strain_naams,
                        2
                   ))
for pairs in cbn:
    if (pairs[0], pairs[1]) in sum_snp_4_2.keys():
        df_mat.loc[pairs[1], pairs[0]] = sum_snp_4_2[(pairs[0], pairs[1])]
    elif (pairs[1], pairs[0]) in sum_snp_4_2.keys():
        df_mat.loc[(pairs[1], pairs[0])] = sum_snp_4_2[(pairs[1], pairs[0])]

df_mat.iloc[1:, 0: -1].to_csv(args.output_csv)
