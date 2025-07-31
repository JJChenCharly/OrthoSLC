from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

import sys#, getopt
import argparse
import os

from OrthoSLC import __version__

argv = sys.argv[1:]

p = argparse.ArgumentParser(
    prog="python3 TK_AlnConcat.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)

p.add_argument("-i", "--input_path",
               metavar='',
                required=True,
                help="<dir> path/to/input/directory")
p.add_argument("-o", "--output_path",
               metavar='',
                required=True,
                help="<dir> path/to/output/file")
p.add_argument("-T", "--ID_TSV",
               metavar='',
                required=True,
                help="<txt> path/to/Step1.txt")
p.add_argument("-f", "--output_format",
               metavar='',
                required=True,
                help="<str> e.g., 'phylip-relaxed', 'fasta'")

args, extras = p.parse_known_args()

with open(args.ID_TSV, 'r') as f_read:
    rls = f_read.readlines()
    
d = {}
for x in rls:
    l = x.split('\t')
    d[l[0]] = l[1]
inv_map = d

to_write = []

seq_coll = {}

# 把fasta 变成字典
ini_seq = SeqIO.to_dict(SeqIO.parse(os.path.join(
     args.input_path, 
     os.listdir(args.input_path)[0]
     ),
     "fasta")
     )
for x in ini_seq:
    minus_ind_in_x = x.index("-")
    
    seq_coll[inv_map[x[0 : minus_ind_in_x]]] = ""

# 
for x in os.listdir(args.input_path):
#     print(x)
    x_fasta = SeqIO.to_dict(SeqIO.parse(os.path.join(args.input_path, x),  
                                        "fasta"))
    for y in x_fasta:
        minus_ind_in_y = y.index("-")
        
        seq_coll[inv_map[y[0 : minus_ind_in_y]
                        ]
                ] = seq_coll[inv_map[y[0 : minus_ind_in_y]
                                    ]
                            ] + str(x_fasta[y].seq)

for z in seq_coll:
    to_write.append(SeqRecord(id = z, 
                              name = z,
                              description = "",
                              dbxrefs=[],
                              seq = Seq(seq_coll[z])
                             )
                   )
with open(args.output_path, "w") as output_handle:
        SeqIO.write(to_write, output_handle, args.output_format)