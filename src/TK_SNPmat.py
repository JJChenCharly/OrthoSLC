from Bio import SeqIO
import pandas as pd

import argparse
from pathlib import Path

from itertools import combinations
from multiprocessing import Pool
from collections import defaultdict

import numpy as np

from OrthoSLC import __version__

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

def _init_worker(base_path, species_order):
    """Store shared configuration for worker processes."""
    global _KALIGNED_PATH, _SPECIES_ORDER, _PAIR_INDEX

    _KALIGNED_PATH = base_path
    _SPECIES_ORDER = species_order
    _PAIR_INDEX = [(i, j) for i in range(len(species_order))
                   for j in range(i + 1, len(species_order))]


def _encode_sequence(sequence: str) -> np.ndarray:
    """Convert a nucleotide sequence into a numpy byte array."""
    return np.frombuffer(sequence.upper().encode("ascii"), dtype="S1")


def _process_alignment(filename: str):
    """Compute pairwise SNP counts for a single alignment file."""
    filepath = Path(_KALIGNED_PATH) / filename
    fasta_records = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))

    sequence_bytes = {}
    for record_id, record in fasta_records.items():
        if "-" not in record_id:
            raise ValueError(
                f"Record id '{record_id}' in '{filepath}' does not contain '-' separator."
            )
        species_id = record_id.split("-", 1)[0]
        sequence_bytes[species_id] = _encode_sequence(str(record.seq))

    missing_species = [sid for sid in _SPECIES_ORDER if sid not in sequence_bytes]
    if missing_species:
        raise KeyError(
            f"Alignment '{filepath}' is missing sequences for: {', '.join(missing_species)}"
        )

    seq_len_set = {len(seq) for seq in sequence_bytes.values()}
    if len(seq_len_set) != 1:
        raise ValueError(
            f"Alignment '{filepath}' contains sequences with inconsistent lengths: {seq_len_set}"
        )

    ordered_matrix = np.vstack([sequence_bytes[sid] for sid in _SPECIES_ORDER])
    diff_matrix = ordered_matrix[:, None, :] != ordered_matrix[None, :, :]
    snp_matrix = diff_matrix.sum(axis=2, dtype=np.int32)

    pair_counts = {
        (_SPECIES_ORDER[i], _SPECIES_ORDER[j]): int(snp_matrix[i, j])
        for i, j in _PAIR_INDEX
    }

    return pair_counts, ordered_matrix.shape[1]

if __name__ == '__main__':
    procc_num = int(args.resource_total)

    kaligned_dir = Path(kaligned_path)
    input_files = sorted(
        [f.name for f in kaligned_dir.iterdir() if f.is_file() and not f.name.startswith('.')]
    )

    if not input_files:
        raise FileNotFoundError(f"No alignment files found under '{kaligned_path}'.")

    first_file = kaligned_dir / input_files[0]
    species_ids = sorted(
        record.id.split("-", 1)[0]
        for record in SeqIO.parse(first_file, "fasta")
    )

    worker_count = max(1, min(procc_num, len(input_files)))

    with Pool(processes=worker_count, initializer=_init_worker,
              initargs=(str(kaligned_dir), species_ids)) as pool:
        results = pool.map(_process_alignment, input_files)

    sum_snp = defaultdict(int)

    for pair_counts, _ in results:
        for key, value in pair_counts.items():
            sum_snp[key] += value

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
                     , dtype=np.int64
                     )
df_mat.index = strain_naams
df_mat.columns = strain_naams

cbn = list(combinations(strain_naams,
                        2
                   ))
for pairs in cbn:
    if (pairs[0], pairs[1]) in sum_snp_4_2.keys():
        value = int(round(sum_snp_4_2[(pairs[0], pairs[1])]))
        df_mat.loc[pairs[1], pairs[0]] = value
        df_mat.loc[pairs[0], pairs[1]] = value
    elif (pairs[1], pairs[0]) in sum_snp_4_2.keys():
        value = int(round(sum_snp_4_2[(pairs[1], pairs[0])]))
        df_mat.loc[pairs[1], pairs[0]] = value
        df_mat.loc[pairs[0], pairs[1]] = value

df_mat.iloc[1:, 0: -1].to_csv(args.output_csv)
